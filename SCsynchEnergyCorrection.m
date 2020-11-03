function [deltaF,ERROR] = SCsynchEnergyCorrection(SC,varargin)
% SCsynchEnergyCorrection
% =======================
%
% NAME
% ----
% SCsynchEnergyCorrection - Calculates a beam based correction to the rf frequency
%
% SYNOPSIS
% --------
% `[deltaF, ERROR] = SCsynchEnergyCorrection(SC [, options])`
%
%
% DESCRIPTION
% -----------
%
% Changes the cavity frequency within the user defined interval and evaluates the mean
% turn-by-turn horizontal BPM deviation. A straight line is fitted to the data and the zero crossing
% is identified. It is assumed that the beam is injected relatively close to the synchronous phase.
% Then, for sufficiently small number of turns the synchrotron motion results in a mean turn-by-turn
% energy deviation which is zero if the synchronous energy defined by by the rf frequency matches
% the injected beam energy. The number of evaluated turns should be smaller than sunchrotron period. 
% Note that if more than one cavity is specified the same frequency steps are
% applied to all cavities. Not also that results might be compromised if the beam transmission
% differs significantly throughout the different applied frequency steps.
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
% `'range'` (`[-1E3,1E3]`)::
%	Frequency range [Hz]
% `'nSteps'` (`15`)::
%	Number of frequency steps to be evaluated
% `'nTurns'` (`150`)::
%	Number of turns to be evaluated
% `'minTurns'` (`50`)::
%	Minimum number of turns the beam has to survive in order to be included in the calculation.
% `'plotResults'` (`0`)::
%	If true, results are plotted.
% `'plotProgress'` (`0`)::
%	If true, each frequency step is plotted.
% `'verbose'` (`0`)::
%	If true, debug information is printed.
%
%
% RETURN VALUES
% -------------
% `deltaF`::
%	Frequency correction step to be added to cavity frequency.
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
	addOptional(p,'range',[-1E3,1E3]);
	addOptional(p,'nSteps',15);
	addOptional(p,'nTurns',150);
	addOptional(p,'minTurns',0);
	addOptional(p,'plotResults',0);
	addOptional(p,'plotProgress',0);
	addOptional(p,'verbose',0);
	parse(p,varargin{:});
	par = p.Results;

	inputCheck(SC,par)

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Initialization

	% Initialize output
	ERROR  = 0;
	deltaF = 0;

	% Define frequency vector
	fTestVec = linspace(par.range(1),par.range(2),par.nSteps);

	% Initialize TbT BPM shift
	BPMshift = nan(1,length(fTestVec));

	% Adjust number of turns
	SC.INJ.nTurns = par.nTurns;

	if par.verbose
		fprintf('Correct energy error with: \n %d Particles \n %d Turns \n %d Shots \n %d Frequency steps between [%.1f %.1f]kHz.\n\n',SC.INJ.nParticles,SC.INJ.nTurns,SC.INJ.nShots,par.nSteps,1E-3*par.range)
	end

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Main script

	% Loop over different frequencies
	for nE = 1:length(fTestVec)

		% Change RF frequency
		tmpSC = SCsetCavs2SetPoints(SC,par.cavOrd,'Frequency',fTestVec(nE),'add');

		% Calculate turn-by-turn horizontal trajectory shift
		[BPMshift(nE),TBTdE] = getTbTEnergyShift(tmpSC,par);

		% Plot current step
		if par.plotProgress
			plotProgress(TBTdE,BPMshift,fTestVec,nE)
		end
	end

	% Define values for fitting
	x = fTestVec(:);
	y = BPMshift(:);

	% Exclude nan
	x(isnan(y)) = [];
	y(isnan(y)) = [];

	if isempty(y)
		ERROR = 1;
		fprintf('No transmission.\n')
		return
	end

	% Fit line
	p = polyfit(x,y,1);

	% Find zero crossing
	deltaF = -p(2)/p(1);

	% Plot results
	if par.plotResults
		plotFunction()
	end

	if isnan(deltaF)
		ERROR = 2;
		fprintf('NaN energy correction step.')
		return
	end


	% Print results
	if par.verbose
		% Calculate reference closed orbit
		XCO = findorbit6(SC.RING);
		% Apply frequency correction
		tmpSC = SCsetCavs2SetPoints(SC,par.cavOrd,'Frequency',deltaF,'add');
		% Calculate final closed orbit
		XCOfinal = findorbit6(tmpSC.RING);
		% Print results
		fprintf('Frequency correction step: %.2fkHz\n',1E-3*deltaF)
		fprintf('>> Energy error corrected from %.2f%% to %.2f%%\n',1E2*(SC.INJ.Z0(5) - XCO(5)),1E2*(SC.INJ.Z0(5) - XCOfinal(5)))
	end

	% Plot final results
	function plotFunction
		figure(88);clf;hold on
		% Plot TBT BPM shift
		plot(1E-3*fTestVec,1E6*BPMshift,'o')
		% Plot fit
		plot(1E-3*fTestVec,1E6*(fTestVec*p(1)+p(2)),'--')
		% Plot correction
		plot(1E-3*(deltaF),0,'kX','MarkerSize',16)
		xlabel('$\Delta f$ [$kHz$]');ylabel('$<\Delta x>$ [$\mu$m/turn]');
		legend({'Measurement','Fit','dE correction'});%,'Closed orbit'})
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
	if par.minTurns>=par.nTurns
		error('Minimum number of turns (%d) must be smaller than overal number of turns (%d).', par.minTurns,par.nTurns)
	end
	if ~isfield(SC.RING{par.cavOrd},'Frequency')
		error('This is not a cavity (ord: %d)',par.cavOrd)
	end
	if ~any(strcmp(SC.RING{par.cavOrd}.PassMethod,{'CavityPass','RFCavityPass'}))
		warning('Cavity (ord: %d) seemed to be switched off.',par.cavOrd)
	end
end

% Calculate mean turn-by-turn horizontal BPM shift
function [BPMshift,TBTdE] = getTbTEnergyShift(SC,par)
	% Calculate beam reading
	B = SCgetBPMreading(SC);

	% Reshape horizontal BPM readings turn-by-turn
	BB = reshape(B(1,:),[],SC.INJ.nTurns);

	% Calculate mean turn-by-turn difference (for plotting)
	TBTdE = mean( BB - repmat(BB(:,1),1,SC.INJ.nTurns) );

	% Calculate mean turn-by-turn energy shift
	x = (1:SC.INJ.nTurns)';
	y = TBTdE';

	% Exclude nan BPM readings
	x(isnan(y)) = [];
	y(isnan(y)) = [];

	% Check  minimum turn requirement
	if length(y)<par.minTurns
		BPMshift = nan;
	else
		% Perform least squares fit
		BPMshift = x\y;
	end
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
	xlabel('$\Delta f$ [Hz]');ylabel('$<\Delta x>$ [m/turn]');
	set(findall(gcf,'-property','FontSize'),'FontSize',18);
	set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
	set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
	set(gcf,'color','w');drawnow
	drawnow
end
