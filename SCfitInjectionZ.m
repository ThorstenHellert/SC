function [deltaZ0,ERROR] = SCfitInjectionZ(SC,mode,varargin)
% SCfitInjectionZ
% ===============
%
% NAME
% ----
% SCfitInjectionZ - Fits the injected beam trajectory
%
% SYNOPSIS
% --------
% `[deltaZ0, ERROR] = SCfitInjectionZ(SC, mode [, options])`
%
%
% DESCRIPTION
% -----------
% This function calculates an transverse injection correction based on the BPM readings. Depending on the specified `mode`
% different approaches are being used. In  
%
%
% INPUTS
% ------
% `SC`::
%	SC base structure.
% `mode`::
%	Method to identify the injection offset. Possible are
%   - `'fitTrajectory'`:: Based on the ideal lattice a trajectory is fitted which best matches the first `N` BPM readings 
%                         as defined in the options.
%   - `'injectionDrift'`:: It is assumed that a 1-turn period orbit is established (see *SCfeedbackBalance*) and that 
%                         between the last BPM of the first turn and the first BPM in the second turn is a drift space.
%                         A linear regression is used to identify the injected beam trajectory. 
%
% OPTIONS
% -------
% The following options can be specified as name-value pairs:
%
% `'nDims'` (1:2)::
%	Which transverse planes should be considered.
% `'nBPMs'` (1:3)::
%	Which BPMs.
% `'plotFlag'` (0)::
%	If true, results are plotted.
% `'verbose'` (0)::
%	If true, additional information is printed.
%
%
% RETURN VALUE
% ------------
% `deltaZ0`::
% 	Injected beam trajectory correction.
% `ERROR`::
% 	Error flag.

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Input check  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% Parse optional arguments
p = inputParser;
addOptional(p,'nDims',1:2);
addParameter(p,'nBPMs',[1:3]);
addOptional(p,'nShots',SC.INJ.nShots);
addOptional(p,'verbose',0);
addOptional(p,'plotFlag',0);
parse(p,varargin{:});
par = p.Results;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Initialization % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Initialize error flag
ERROR = 0;

% Initialize output
deltaZ0 = zeros(6,1);

inputCheck()
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Main function  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Adjust number of shots
SC.INJ.nShots = par.nShots;

% Get BPM reading
B = SCgetBPMreading(SC);

switch mode
	case 'fitTrajectory'

		% Get BPM ordinates
		ordsUsed = SC.ORD.BPM(par.nBPMs);

		% Get reference BPM reading
		Bref = B(:,par.nBPMs);

		% Search for optimal injection trajectory
		deltaZ0(1:4) = -fminsearch(@merritFunction,zeros(4,1));

		% Plot results
		if par.plotFlag
			% Aplly correction
			SC.INJ.Z0 = SC.INJ.Z0 + deltaZ0;
			B1 = SCgetBPMreading(SC);
			% Get BPM s-Positions
			sBPM = findspos(SC.RING,SC.ORD.BPM(par.nBPMs));
			figure(342);clf;titleStr={'Horizontal','Vertical'};
			for nDim=1:2
				subplot(1,2,nDim)
				plot(sBPM,1E3*Bref(nDim,:),'O--',sBPM,1E3*B1(nDim,par.nBPMs),'X--')
				xlabel('s [m]');ylabel('BPM reading [mm]');title(titleStr{nDim});legend({'Initial','After correction'})
			end
			set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
			set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
			set(findall(gcf,'-property','FontSize'),'FontSize',16);
			set(gcf,'color','w');
			drawnow
		end
		
	case 'injectionDrift'

		% Get BPM s-Positions
		tmpS = findspos(SC.RING,SC.ORD.BPM);
		sBPM = [tmpS(end) - findspos(SC.RING,length(SC.RING)+1) , tmpS(1)];

		% Get BPM readings at injection section between 1st and 2nd turn
		Bref = [B(:,length(SC.ORD.BPM)) B(:,length(SC.ORD.BPM)+1)];

		% Loop over planes
		for nDim=par.nDims
			% Perform line fit
			sol{nDim} = polyfit(sBPM,Bref(nDim,:),1);

			% Write injection correction
			deltaZ0(2*nDim-1) = - sol{nDim}(2);
			deltaZ0(2*nDim)   = - sol{nDim}(1);
		end
		
		% Plot results
		if par.plotFlag
			figure(342);clf;titleStr={'Horizontal','Vertical'};
			for nDim=par.nDims
				subplot(1,2,nDim)
				plot(sBPM,1E6*Bref(nDim,:),'o',sBPM,1E6*(sol{nDim}(1)*sBPM+sol{nDim}(2)),'--',sBPM,1E6*(SC.INJ.Z0(2*nDim)*sBPM+SC.INJ.Z0(2*nDim-1)),'k-',sBPM,[0 0],'k--')
				legend({'BPM reading','Fitted trajectory','Real trajectory'});xlabel('s [m]');ylabel('Beam offset [mm]');title(titleStr{nDim});
			end
			set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
			set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
			set(findall(gcf,'-property','FontSize'),'FontSize',16);
			set(gcf,'color','w');
			drawnow
		end

end

% Print results
if par.verbose
	fprintf('\nInjection trajectory corrected from \n x:  %.0fum -> %.0fum \n x'': %.0furad -> %.0furad \n y:  %.0fum -> %.0fum \n y'': %.0furad -> %.0furad\n',...
		1E6*SC.INJ.Z0(1),1E6*(SC.INJ.Z0(1)+deltaZ0(1)),1E6*SC.INJ.Z0(2),1E6*(SC.INJ.Z0(2)+deltaZ0(2)),1E6*SC.INJ.Z0(3),1E6*(SC.INJ.Z0(3)+deltaZ0(3)),1E6*SC.INJ.Z0(4),1E6*(SC.INJ.Z0(4)+deltaZ0(4)))
end

% Check for errors
if any(isnan(deltaZ0))
	ERROR = 1;
end

function inputCheck()
	switch mode
		case 'injectionDrift'
			if SC.INJ.nTurns~=2
				error('Injection pattern (''SC.INJ.nTurns'') must be two turns.')
			end
		case 'fitTrajectory'
			if min(diff(par.nBPMs))<=0
				error('BPMs must be given in correct order.')
			end
		otherwise
			error('Unsupported mode: ''%s''',mode)
	end
end

% Merrit function for tracking
function out = merritFunction(x)
	% Get trajectory from ideal lattice
	Ta = atpass(SC.IDEALRING, [x;0;0], 1, 1, ordsUsed);
	% Copy horizontal and vertical offsets
	T = Ta([1 3],:);
	% Calculate rms value
	out = sqrt(mean((Bref(:)-T(:)).^2));
end
	
	
end
