function [deltaPhi,ERROR] = SCsynchPhaseCorrection(SC,varargin)
% SCsynchPhaseCorrection
% ======================
%
% NAME
% ----
% SCsynchPhaseCorrection - Calculates a beam based correction to the rf cavity phase
%
% SYNOPSIS
% --------
% `[deltaPhi, ERROR] = SCsynchPhaseCorrection(SC [, options])`
%
%
% DESCRIPTION
% -----------
%
% Changes the cavity phase within the phase interval [-pi,pi] stepwise and evaluates the mean
% turn-by-turn horizontal BPM deviation. A sinusoidal function is fitted to the data and the zero
% crossing is identified. It is assumed that for sufficiently small number of turns the synchrotron
% motion is negligible and injection at the synchronous phase will result in zero turn-by-turn
% energy variation. Thus, the horizontal turn-by-turn BPM readings should in first approximation not
% differ. The number of evaluated turns should be significantly smaller than half a sunchrotron
% period. Note that if more than one cavity is specified the same phase steps are applied to all
% cavities. Results might be compromised.
%
%
% INPUTS
% ------
% `SC`::
%   The SC base structure
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'cavOrd'` (`SC.ORD.Cavity`)::
%	Ordinate of evaluated cavity
% `'nSteps'` (`30`)::
%	Number of phase steps to be evaluated
% `'nTurns'` (`30`)::
%	Number of turns to be evaluated
% `'plotResults'` (`0`)::
%	If true, results are plotted.
% `'plotProgress'` (`0`)::
%	If true, each phase step is plotted.
% `'verbose'` (`0`)::
%	If true, debug information is printed.
%
%
% RETURN VALUES
% -------------
% `deltaPhi`::
%	Phase correction step to be added to cavity field 'TimeLag'.
% `ERROR`::
%	Error value.
%
% ERRORS
% ------
% `0`::
% 	All good.
% `1`::
%	Horizontal TBT BPM deviation shows no zero crossing
% `2`::
% 	Sinusoidal fit function shows no zero crossing
%
%
% SEE ALSO
% --------
% *SCsetCavs2SetPoints*, *SCgetBPMreading*, *SCsynchEnergyCorrection*


	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Input check

	% Parse input
	p = inputParser;
	addOptional(p,'cavOrd',SC.ORD.Cavity);
	addOptional(p,'nSteps',15);
	addOptional(p,'nTurns',20);
	addOptional(p,'plotResults',0);
	addOptional(p,'plotProgress',0);
	addOptional(p,'verbose',0);
	parse(p,varargin{:});
	par = p.Results;

	inputCheck(SC,par);

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Initialization

	% Initialize output
	ERROR    = 0;
	deltaPhi = 0;

	% Initialize TbT BPM shift
	BPMshift = nan(1,par.nSteps);

	% RF wavelength [m]
	lambda = 299792458/SC.RING{par.cavOrd}.Frequency;

	% Define phase scan range
	lambdaTestVec = 1/2 * lambda * linspace(-1,1,par.nSteps);

	% Adjust number of turns
	SC.INJ.nTurns = par.nTurns;

	if par.verbose
		fprintf('Calibrate RF phase with: \n %d Particles \n %d Turns \n %d Shots \n %d Phase steps.\n\n',SC.INJ.nParticles,SC.INJ.nTurns,SC.INJ.nShots,par.nSteps)
	end
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Main script


	% Loop over different phase offsets
	for nL = 1:length(lambdaTestVec)

		% Change phase
		tmpSC = SCsetCavs2SetPoints(SC,par.cavOrd,'TimeLag',lambdaTestVec(nL),'add');
	
		% Calculate turn-by-turn horizontal trajectory shift
		[BPMshift(nL),TBTdE] = getTbTEnergyShift(tmpSC);
		
		
		% Plot current step
		if par.plotProgress
			plotProgress(TBTdE,BPMshift,lambdaTestVec,nL)
		end
	end

	% Check if zero crossing is reached in data
	if ~(max(BPMshift)>0 && min(BPMshift)<0)
		fprintf('Zero crossing not within data set.\n')
		ERROR = 1;
		return
	end

	% Define sinusoidal function
	sinFun = @(par,s) par(1)*sin(2*pi*(par(4)*s + par(2))) +par(3);

	% Define merit function
	fomFun = @(par)  sum( (sinFun(par,lambdaTestVec) - BPMshift).^2 );

	% Loop over different start point guesses
	for startPhase=[-pi,-pi/2,-pi/4,0]
		% Run fminsearch to determine sinusoidal parameters
		sol = fminsearch(fomFun, [max(BPMshift)-mean(BPMshift);  startPhase;  mean(BPMshift) ; 0.5/max(lambdaTestVec)]);

		% Save solution in sinusoidal function
		solFun = @(x) sinFun(sol,x);

		% Generate densely populated lambda-vector
		xp = linspace(lambdaTestVec(1),lambdaTestVec(end),100);

		% Check if zero crossing is reached in fitted function
		if ~(max(solFun(xp))>0 && min(solFun(xp))<0)
			fprintf('Zero crossing not within fit function, trying different start point guess.\n')
		else
			break
		end
	end

	% Check if zero crossing is reached in fitted function
	if ~(max(solFun(xp))>0 && min(solFun(xp))<0)
		fprintf('Zero crossing not within fit function\n')
		ERROR = 2;
% 		return
	end

	% Identify maximum value
	[~,maxValInd] = max(solFun(xp));

	% Find zero crossing on the left hand side of maximum value
	deltaPhi = fzero(solFun,xp(maxValInd)-abs(lambdaTestVec(1))/2);

	% Shift phase into [-pi,pi]-intervall
	if deltaPhi > lambda/2
		deltaPhi = deltaPhi - lambda;
	elseif deltaPhi < -lambda/2
		deltaPhi = deltaPhi + lambda;
	end

	% Plot results
	if par.plotResults
		plotFunction()
	end

	% Check for NaN results
	if isnan(deltaPhi)
		ERROR = 3;
		fprintf('SCrfCommissioning: ERROR (NaN phase)\n')
		return
	end

	% Print results
	if par.verbose
		% Initial phase
		XCO   = findorbit6(SC.RING);
		% Apply phase correction
		tmpSC = SCsetCavs2SetPoints(SC,par.cavOrd,'TimeLag',deltaPhi,'add');
		% Calculate corrected closed orbit
		XCOfinal = findorbit6(tmpSC.RING);
		% Calculate figures of merrit
		initial = rem((XCO(6)     -SC.INJ.Z0(6))/lambda*360,360);
		final   = rem((XCOfinal(6)-SC.INJ.Z0(6))/lambda*360,360);
		% Print result
		fprintf('Phase correction step: %.3m\n',deltaPhi)
		fprintf('>> Static phase error corrected from %.0fdeg to %.1fdeg\n',initial,final)
	end

	% Plot final results
	function plotFunction
		figure(87);clf;hold on
		% Plot TBT BPM shift
		plot((lambdaTestVec+SC.INJ.Z0(6))/lambda*360,1E6*BPMshift,'o')
		% Plot fitted sinosoidal functions
		plot((xp+SC.INJ.Z0(6))/lambda*360,1E6* solFun(xp))
		% Plot initial phase
		plot((SC.INJ.Z0(6))/lambda*2*180,1E6* (sol(1)*sin(2*pi*(sol(4)*0 + sol(2))) +sol(3)),'rD','MarkerSize',16)
		% Plot final phase
		plot((deltaPhi)/lambda*2*180,0,'kX','MarkerSize',16)
		set(gca,'XLim',180*[-1 1],'box','on');
		legend({'Measurement','Fit','Initial phase','Phase correction'})
		xlabel('RF phase [$^\circ$]');ylabel('$<\Delta x>$ [$\mu$m/turn]');
		set(findall(gcf,'-property','FontSize'),'FontSize',18);
		set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
		set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
		set(gcf,'color','w');drawnow
	end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions

% Input check
function inputCheck(SC,par)
	if ~strcmp(SC.INJ.trackMode,'TBT')
		error('Trackmode should be turn-by-turn (''SC.INJ.trackmode=TBT'').')
	end
	if par.nSteps<2 || par.nTurns<2
		error('Number of steps and number of turns must be larger than 2.')
	end
	if ~isfield(SC.RING{par.cavOrd},'Frequency')
		error('This is not a cavity (ord: %d)',par.cavOrd)
	end
	if ~any(strcmp(SC.RING{par.cavOrd}.PassMethod,{'CavityPass','RFCavityPass'}))
		warning('Cavity (ord: %d) seemed to be switched off.',par.cavOrd)
	end
end
% Calculate mean turn-by-turn horizontal BPM shift
function [BPMshift,dE] = getTbTEnergyShift(SC)
	% Calculate beam reading
	B = SCgetBPMreading(SC);

	% Reshape horizontal BPM readings turn-by-turn
	BB = reshape(B(1,:),[],SC.INJ.nTurns);

	% Calculate mean turn-by-turn difference
	dE = mean( BB - repmat(BB(:,1),1,SC.INJ.nTurns) );

	% Calculate mean turn-by-turn energy shift
	x = (1:(SC.INJ.nTurns-1))';
	y = dE(2:end)';

	% Exclude nan BPM readings
	x(isnan(y),:) = [];
	y(isnan(y))   = [];

	% Perform least squares fit
	BPMshift= x\y;
end

% Plot current step
function plotProgress(TBTdE,BPMshift,fTestVec,nE)
	figure(2);clf
	subplot(2,1,1);hold on
	plot(TBTdE,'o');
	plot([1:length(TBTdE)] * BPMshift(nE),'--')
	xlabel('Number of turns');ylabel('$<\Delta x_\mathrm{TBT}>$ [m]');
	subplot(2,1,2);
	plot(fTestVec(1:nE),BPMshift(1:nE),'o')
	xlabel('$\Delta \phi$ [m]');ylabel('$<\Delta x>$ [m/turn]');
	set(findall(gcf,'-property','FontSize'),'FontSize',18);
	set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
	set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
	set(gcf,'color','w');drawnow
	drawnow
end
