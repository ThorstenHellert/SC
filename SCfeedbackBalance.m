function [SC,ERROR] = SCfeedbackBalance(SC,Mplus,varargin)
% SCfeedbackBalance
% =================
%
% NAME
% ----
% SCfeedbackBalance - balance two-turn BPM readings
%
% SYNOPSIS
% --------
% `[SC, ERROR] = SCfeedbackBalance(SC, Mplus [,options])`
%
%
% DESCRIPTION
% -----------
% Generates a period-1 closed orbit, after two-turn transmission has been
% achieved. This is done by iteratively applying correction steps, calculated
% based on the pseudo-inverse two-turn trajectory response matrix `Mplus`.  The
% trajectory in the first turn is corrected towards the reference orbit `R0`,
% whereas the trajectory in the second turn is corrected towards the trajectory
% measured in the first turn; this approach seems to be more stable than the
% directly application of the two-turn TRM to the two-turn BPM readings.
% It converges to a state where BPM readings in both turns are very similar,
% indicating a period-1 closed orbit.
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
% `'eps'` (`1e-5`)::
%	break, if the coefficient of variation of the RMS BPM reading is below
%	this value
% `'R0'` (zeros)::
%	target orbit in the format `[x_1 ... x_n y_1 ...y_n]`, where
%	`[x_i,y_i]` is the target position at the i-th BPM.
% `'maxsteps'` (10)::
%	A maximum of `'maxsteps'` correction steps is performed.
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
%
% ERRORS
% ------
% `0`::
% 	All fine.
% `1`::
%	A correction step resulted in less transmission, than before.
% `2`::
% 	Transmission was lost during correction procedure. This is an indicator
% 	that `Mplus` might be unstable.
% `3`::
%	The feedback was unstable, when 'maxsteps' was reached.
%
% SEE ALSO
% --------
%  *SCgetPinv*, *SCfeedbackRun*, *SCfeedbackFirstTurn*, *SCfeedbackStitch*, *SCgetBPMreading*, *SCsetCMs2SetPoints*


	% Parse optional arguments
	p = inputParser;
	addOptional(p,'eps',1E-5);
	addOptional(p,'R0',zeros(2,1));
	addOptional(p,'maxsteps',10);
	addOptional(p,'CMords',SC.ORD.CM);
	addOptional(p,'BPMords',SC.ORD.BPM);
	addOptional(p,'verbose',0);
	parse(p,varargin{:});
	par=p.Results;


	if par.verbose; fprintf('SCfeedbackBalance: Start\n'); end;

	% Initialize error state
	ERROR = 0;

	% Initialize history. ``hist'' stores the 100 last last-reached BPM indices.
	% Based on this history, some break conditions are checked.
	hist = -1*ones(1,100);

	% Initialize history. ``BRMShist'' stores the 100 last RMS BPM readings.
	% Based on this history, some break conditions are checked.
	BRMShist = nan(1,100);

	% Main loop. In each step, a correction step is applied and the
	% break-conditions are checked.
	cnt=1;
	for steps=1:par.maxsteps

		B = SCgetBPMreading(SC,'BPMords',par.BPMords); % Inject ...
		correctionStep(); % call correction subroutine. See below.

		if isSetback(hist)
			% Fail, is we have less transmission than before.
			if par.verbose; fprintf('SCfeedbackBalance: FAIL (setback)\n'); end;
			ERROR = 1; return;
		end

		if ~isTransmit(hist)
			% Fail, if we do not have full transmission.
			if par.verbose; fprintf('SCfeedbackBalance: FAIL (lost transmission)\n'); end;
			ERROR = 2; return;
		else
			ERROR = 0;
		end

		if isConverged(BRMShist,3,par.eps)
			% Succeed, if the RMS BPM-reading has converged.
			if par.verbose; fprintf('SCfeedbackBalance: Success (converged after %d steps)\n',steps); end;
			ERROR = 0; return;
		end
		cnt = cnt+1;
	end

	% This part is reached only if more than 'maxsteps' steps have been
	% performed.

	if isStable(min(10,par.maxsteps),par.eps)
		% Reaching 'maxsteps' with the feedback being stable is
		% considered a success.
		if par.verbose; fprintf('SCfeedbackBalance: Success (maxsteps reached)\n'); end;
		ERROR=0; return
	else
		% Reaching 'maxsteps' with the feedback being unstable is
		% considered a failure.
		if par.verbose; fprintf('SCfeedbackBalance: FAIL (maxsteps reached, unstable)\n'); end;
		ERROR=3; return
	end


%%%%%%%%%%%%%%%% SUBROUTINES

	function correctionStep()
			% Sub routine that calculates the correction step by
			% applying the pseudo-inverse matrix to the current
			% first-turn orbit deviation from the target-orbit 'R0'
			% and second-turn deviation from the first turn.
			% The correction step is then applied and the old
			% RMS-beam reading is prepended to 'BRMShist'. The last
			% reached BPM index is stored in hist.
			hist = logLastBPM(hist,B);
			lBPM = size(B,2);
			Bx1 = B(1,1:lBPM/2);
			By1 = B(2,1:lBPM/2);
			Bx2 = B(1,(lBPM/2+1):end);
			By2 = B(2,(lBPM/2+1):end);
			DELTABx=Bx2-Bx1;
			DELTABy=By2-By1;
			R = [Bx1-par.R0(1,:),DELTABx,By1-par.R0(2,:),DELTABy]';
			R(isnan(R))=0;
			dphi = Mplus * R;
			BRMShist = circshift(BRMShist,1);
			BRMShist(1) = sqrt(var(R,1));
			SC = SCsetCMs2SetPoints(SC,par.CMords{1},-dphi(1:length(par.CMords{1})      ),1,'add');
			SC = SCsetCMs2SetPoints(SC,par.CMords{2},-dphi(  length(par.CMords{1})+1:end),2,'add');
	end



	function res = isStable(n,eps)
		% Is true, if the coefficient of variation of the last 'n' RMS
		% BPM-readings is less that 'eps'.
		CV=var(BRMShist(1:n),1)/std(BRMShist(1:n));
		res = CV<eps;
	end

end





%%%%%%%%%%%%%%%% Helper Functions

function res = isSetback(hist)
	% True, if transmission is less than before
	res = hist(1)<hist(2);
end

function res = isTransmit(hist)
	% True, if we have full transmission
	res = hist(1)==0;
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

function res = isConverged(hist,n,eps)
	% Preliminary convergence criterion
	% TODO: This is not a very good convergence criterion. PhA
	CV=var(hist(1:n),1)/std(hist(1:n));
	res = CV<eps;
end
