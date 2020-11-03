function [SC,ERROR] = SCfeedbackFirstTurn(SC,Mplus,varargin)
% SCfeedbackFirstTurn
% ===================
%
% NAME
% ----
% SCfeedbackFirstTurn - achieves one-turn transmission
%
% SYNOPSIS
% --------
% `[SC, ERROR] = SCfeedbackFirstTurn(SC, Mplus [, options])`
%
%
% DESCRIPTION
% -----------
% Achieves a first turn in `SC.RING`.  This algorithm is based on the idea that
% repeated orbit corrections calculated via a suitably regularized
% pseudo-inverse trajectory-response matrix `Mplus` will drive the BPM readings
% and CM settings to a fixed point.
%
% lim      B_n = const. , with B    = Phi(Mplus . B ),     (1)
%    n->oo                      n+1                n
%
% where mapping Phi maps corrector kicks to BPM-readings.
% The RMS-values of both, BPM readings and CM settings, are determined by the
% regularization of Mplus.  Successively - during the course of repeated
% application of (1) - more and more transmission is achieved throughout the
% ring, more magnets are traversed near their magnetic center (wich is hopefuly
% at least somewhere near the BPM zero-point), resulting in decreased kicks.
% If, however, the beam encounters a heavily displaced quadrupole magnet this
% approach is bound to fail as correction towards the center of the last
% reached BPM does no good, really. In this case the magnet has to be cleared
% using other means than linear algebra.  In this approach the kicks of an
% increasing number of the last reached CMs are deterministialy ``wiggled''
% until transmission to the next BPM is achieved. Then, application of (1) is
% resumed.
%
%
% INPUT
% -----
% `SC`::
%	SC base structure.
% `Mplus`::
%	Pseudo-inverse trajectory-response matrix.
%
% OPTIONS
% -------
% The following options can be specified as name-value pairs:
%
% `'maxsteps'` (`100`)::
%	break, if this number of correction steps have been performed
% `'R0'` (zeros)::
%	target orbit in the format `[x_1 ... x_n y_1 ...y_n]`, where
%	`[x_i,y_i]` is the target position at the i-th BPM.
% `'wiggleAfter'` (`20`)::
%	Number of iterations wihtout increased transmission to start wiggling.
% `'wiggleSteps'` (`64`)::
%	Number of wiggle steps to perform, before incresing the number
%	of wiggling-CMs.
% `'wiggleRange'` (`[50E-6 200E-6]`)::
%	Range ([min,max] in rad) within which to wiggle the CMs.
% `'CMords'` (`SC.ORD.CM`):: 
%   List of CM ordinates to be used for correction
% `'BPMords'` (`SC.ORD.BPM`):: 
%   List of BPM ordinates at which the reading should be evaluated
% `'verbose'` (0)::
%	If true, debug information is printed.
%
%
% RETURN VALUES
% -------------
% `SC`::
%	SC-structure with corrected `SC.RING`.
% `ERROR`::
%	Error value.
%
% ERRORS
% ------
% `0`::
% 	All fine.
% `1`::
%	'maxsteps' was reached, without producing full transmission.
%
% SEE ALSO
% --------
%  *SCgetPinv*, *SCfeedbackRun*, *SCfeedbackStitch*, *SCfeedbackBalance*, *SCgetBPMreading*, *SCsetCMs2SetPoints*


	% Parse optional arguments
	p = inputParser;
	addOptional(p,'R0',zeros(size(Mplus,2),1));
	addOptional(p,'maxsteps',100);
	addOptional(p,'wiggleAfter',20);
	addOptional(p,'wiggleSteps',64);
	addOptional(p,'wiggleRange',[50E-6 200E-6]);
	addOptional(p,'CMords',SC.ORD.CM);
	addOptional(p,'BPMords',SC.ORD.BPM);
	addOptional(p,'verbose',0);
	parse(p,varargin{:});
	par=p.Results;

	if par.verbose; fprintf('SCfeedbackFirstTurn: Start\n'); end;

	%Initialize history, error value and number of wiggled CMs.
	hist = -1*ones(1,100);
	ERROR = 1;
	nWiggleCM=1;

	% Main loop. In each step, a correction step is applied and the
	% break-conditions are checked.
	for i=1:par.maxsteps

		B = SCgetBPMreading(SC,'BPMords',par.BPMords); % Inject...
		correctionStep(); % call correction subroutine.

		if isRepro(hist,5) && isTransmit(hist)
			% If we had full transmission in the last 5 injections
			% we succeed.
			if par.verbose; fprintf('SCfeedbackFirstTurn: Success\n'); end;
			ERROR = 0; return;
		elseif isRepro(hist,par.wiggleAfter)
			% If the last-reached BPM has not changed over 20 injections,
			% we start the wiggling procedure.
			if par.verbose; fprintf('SCfeedbackFirstTurn: Wiggling\n'); end;

			% Determine last-reached BPMs
			CMidxsH = getLastCMsDim(par,B,1,nWiggleCM); % Last CMs in horz
			CMidxsV = getLastCMsDim(par,B,2,nWiggleCM); % Last CMs in vert
			CMordsH = par.CMords{1}(CMidxsH);
			CMordsV = par.CMords{2}(CMidxsV);

			% Calculate wiggling steps, based on a golden-ratio
			% based covering of the 2D-torus.
			pts = [[0;0] goldenDonut(par.wiggleRange(1),par.wiggleRange(2),par.wiggleSteps)];
			dpts = diff(pts,1,2);

			% Loop over wiggling steps.
			for i=1:size(dpts,2)
				% Add wiggling step.
				SPH = dpts(1,i) * ones(length(CMordsH),1); % Horizontal setpoint change
				SPV = dpts(2,i) * ones(length(CMordsV),1); % Vertical setpoint change
				SC = SCsetCMs2SetPoints(SC,CMordsH,SPH,1,'add');
				SC = SCsetCMs2SetPoints(SC,CMordsV,SPV,2,'add');

				% Inject
				W = SCgetBPMreading(SC,'BPMords',par.BPMords);
				% See if we have reached a new BPM
				hist = logLastBPM(hist,W);
				if isNew(hist)
					% Inject 3 additional times
					hist = logLastBPM(hist,SCgetBPMreading(SC,'BPMords',par.BPMords));
					hist = logLastBPM(hist,SCgetBPMreading(SC,'BPMords',par.BPMords));
					hist = logLastBPM(hist,SCgetBPMreading(SC,'BPMords',par.BPMords));
					% Check if all 3 injections had transmission
					% to the newly reached BPM.
					if isRepro(hist,3)
						hist(1:3) = -1; % void last hist
						nWiggleCM = 0; % Reset Wiggler CM number
						break % Continue with feedback
					end
				end
			end

			% This part is reached only if 'wiggleSteps' have been performed,
			% without reaching a new BPM.

			% Apply three feedback steps and then continue wiggling with an
			% increased number of wiggled CMs.
			B = SCgetBPMreading(SC,'BPMords',par.BPMords);
			correctionStep();
			B = SCgetBPMreading(SC,'BPMords',par.BPMords);
			correctionStep();
			B = SCgetBPMreading(SC,'BPMords',par.BPMords);
			correctionStep();
			nWiggleCM=nWiggleCM+1;

		end
	end

	% If we have reached 'maxsteps' without producing full transmission,
	% that is considered a failure.
	if par.verbose; fprintf('SCfeedbackFirstTurn: FAIL (maxsteps reached)\n'); end;
	ERROR = 1; return;



%%%%%%%%%%%%%%%% SUBROUTINES

	function correctionStep()
		% Sub routine that calculates the correction step by applying
		% the pseudo-inverse matrix to the current orbit deviation from
		% the target-orbit 'R0'. The correction step is then applied
		% and the last reached BPM is is prepended to 'hist'.
		hist = logLastBPM(hist,B);
		R = [B(1,:)'; B(2,:)'];
		dR = R-par.R0;
		dR(isnan(dR)) = 0;
		dphi = Mplus * (dR);
		lastCMh = getLastCMsDim(par,B,1,1);
		lastCMv = getLastCMsDim(par,B,2,1);
		dphi(lastCMh+1:length(par.CMords{1})) = 0;
		dphi((length(par.CMords{1})+1)+lastCMv:end) = 0;
		SC = SCsetCMs2SetPoints(SC,par.CMords{1},-dphi(1:length(par.CMords{1})      ),1,'add');
		SC = SCsetCMs2SetPoints(SC,par.CMords{2},-dphi((length(par.CMords{1})+1):end),2,'add');
	end

end

function out = goldenDonut(r0, r1, Npts)
	% Generates Npts points homogeneously filling the annulus with radii r0 and r1.
	% Points are produced deterministically with their radii going monotonically from r0 to r1.
	% input: int r0 (Start radius)
	% input: int r1 (End radius)
	% input: int Npts (Number of points)
	% returns: 2xNpts double array
	out = zeros(2,Npts); % initialize output array
	phi = 2*pi/((1+sqrt(5))/2); % golden ratio
	theta = 0;
	for n=0:Npts-1
		out(:,n+1) = sqrt((r1^2-r0^2)*n/(Npts-1) + r0^2) * [cos(theta); sin(theta)];
		theta = theta + phi;
	end
end

function idxs = getLastCMsDim(par,B,dim,n)
	% Returns the last reached BPM, in dimension dim
	lastBPMidx = find(~isnan(B(1,:)),1,'last');
	if isempty(lastBPMidx)
		lastBPMidx = length(par.BPMords); % the last one
	end
	lastBPMord = par.BPMords(lastBPMidx);
	idxs = find(par.CMords{dim} <= lastBPMord,n,'last');
end

function hist = logLastBPM(hist,B)
	% Write last reached BPM to 'hist'.
	% '0' indicates all BPMs have been reached.
	hist = circshift(hist,1);
	ord = getLastBPMord(B);
	if ord
		hist(1)=ord;
	else
		hist(1)=0;
	end
end

function ord = getLastBPMord(B)
	% Returns the last reached BPM
	ord = find(isnan(B),1)-1;
end

function res = isRepro(hist,N)
	% True, if last-reached BPM is the same for N injections
	res = all(hist(1:N)==hist(1));
end

function res = isTransmit(hist)
	% True, if we have full transmission
	res = hist(1)==0;
end

function res = isNew(hist)
	res = hist(1)~=hist(2);
end
