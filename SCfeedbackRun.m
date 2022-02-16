function [SC,ERROR] = SCfeedbackRun(SC,Mplus,varargin)
% SCfeedbackRun
% =============
%
% NAME
% ----
% SCfeedbackRun - iterative orbit correction
%
% SYNOPSIS
% --------
% `[SC, ERROR] = SCfeedbackRun(SC, Mplus [, options])`
%
%
% DESCRIPTION
% -----------
% Iteratively applies orbit corrections using the pseudoinverse of the
% trajcectory response matrix `Mplus`, until a break-condition specified by one
% of the 'OPTIONS' is met.  
% The dispersion can be included, thus the rf frequency as a correction 
% parameter. If the dispersion is to be included, `Mplus` has to have the size 
% `(length(SC.ORD.CM{1}) + length(SC.ORD.CM{2}) + 1) x length(SC.ORD.BPM)`, otherwise the size 
% `(length(SC.ORD.CM{1}) + length(SC.ORD.CM{2})) x length(SC.ORD.BPM)`, or correspondingly if the CM
% and/or BPM ordinates for the correction is explicitly given (see options below). `SC.RING` is 
% assumed to be a lattice with transmission through all considered turns. This routine will
% also return, if transmssion is lost.
%
%
% INPUT
% -----
% `SC`::
%	SC base structure.
% `Mplus`::
%	Pseudo-inverse trajectory/orbit response matrix.
%
% OPTIONS
% -------
% The following options can be specified as name-value pairs:
%
% `'target'` (`0`)::
%	break, if the RMS BPM reading reaches this value
% `'maxsteps'` (`30`)::
%	break, if this number of correction steps have been performed
% `'eps'` (`1e-5`)::
%	break, if the coefficient of variation of the RMS BPM reading is below
%	this value
% `'R0'` (zeros(size(Mplus,2),1))::
%	target orbit in the format `[x_1 ... x_n y_1 ...y_n]`, where
%	`[x_i,y_i]` is the target position at the i-th BPM.
% `'scaleDisp'` (0)::
%   Scaling factor for and flag indicating if the dispersion is included in the respoinse matrix
% `'CMords'` (`SC.ORD.CM`):: 
%   List of CM ordinates to be used for correction
% `'BPMords'` (`SC.ORD.BPM`):: 
%   List of BPM ordinates at which the reading should be evaluated
% `'weight'` (ones(size(Mplus,2),1))::
%	weighting vector to be used in the CM step calculation in the format `[x_1 ... x_n y_1 ...y_n]`, where
%	`[x_i,y_i]` is the weight at the i-th BPM.
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
%	The feedback was unstable, when 'maxsteps' was reached.
% `2`::
% 	Transmission was lost during correction procedure. This is an indicator
% 	that `Mplus` might be unstable.
%
% EXAMPLES
% --------
%
% Switch to orbit mode, get the model response matrix and dispersion. Calculate the psudo-inverse 
% while scaling the dispersion by 1E7 and using a Tikhonov regularization parameter of 10. 
% Finally, apply  and apply orbit feedback including dispersion.
% ------------------------------------------------------------------
% SC.INJ.trackMode = 'ORB';
% MCO = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'trackMode','ORB');
% eta = SCgetModelDispersion(SC,SC.ORD.BPM,SC.ORD.Cavity);
% MinvCO = SCgetPinv([MCO 1E7*eta],'alpha',10);
% [CUR,ERROR] = SCfeedbackRun(SC,MinvCO,'scaleDisp',1E7);
% ------------------------------------------------------------------
%
% SEE ALSO
% --------
%  *SCgetPinv*, *SCfeedbackStitch*, *SCfeedbackFirstTurn*, *SCfeedbackBalance*, *SCgetBPMreading*, *SCsetCMs2SetPoints*


	% Parse optional arguments
	p = inputParser;
	addOptional(p,'R0',zeros(size(Mplus,2),1));
	addOptional(p,'eps',1e-5);
	addOptional(p,'target',0);
	addOptional(p,'maxsteps',30);
	addOptional(p,'scaleDisp',0);
	addOptional(p,'CMords',SC.ORD.CM);
	addOptional(p,'BPMords',SC.ORD.BPM);
	addOptional(p,'weight',ones(size(Mplus,2),1));
	addOptional(p,'verbose',0);
	parse(p,varargin{:});
	par=p.Results;
	

	if par.verbose; fprintf('SCfeedbackRun: Start\n'); end;

	% Initialize error state
	ERROR = 1;

	% Initialize history. ``hist'' stores the 100 last RMS BPM readings.
	% Based on this history, the break conditions are checked.
	BPMhist = nan(1,100);

	% Main loop. In each step, a correction step is applied and the
	% break-conditions are checked.
	cnt=1;
	for steps=1:par.maxsteps

		B = SCgetBPMreading(SC,'BPMords',par.BPMords); % Inject ...
		correctionStep(); ... call correction subroutine. See below.

		if any(isnan(B(1,:)))
			% Fail, if we lose transmission.
			if par.verbose; fprintf('SCfeedbackRun: FAIL (lost transmission)\n'); end;
			ERROR = 2; return;
		end

		if BPMhist(1)<par.target && isStable(min(10,par.maxsteps),par.eps)
			% Succeed, if we are below our RMS BPM-reading 'target'
			% and the feedback has been stable for the last 10
			% injections.
			if par.verbose; fprintf('SCfeedbackRun: Success (target reached)\n'); end;
			ERROR = 0; return;
		end

		if isConverged(3,par.eps)
			% Succeed, if the RMS BPM-reading has converged.
			if par.verbose; fprintf('SCfeedbackRun: Success (converged after %d steps)\n',steps); end;
			ERROR = 0; return;
		end

		cnt = cnt+1;
	end


	% This part is reached only if more than 'maxsteps' steps have been
	% performed.

	if isStable(min(10,par.maxsteps),par.eps) || par.maxsteps==1
		% Reaching 'maxsteps' with the feedback being stable is
		% considered a success.
		if par.verbose; fprintf('SCfeedbackRun: Success (maxsteps reached)\n'); end;
		ERROR=0; return
	else
		% Reaching 'maxsteps' with the feedback being unstable is
		% considered a failure.
		if par.verbose; fprintf('SCfeedbackRun: FAIL (maxsteps reached, unstable)\n'); end;
		ERROR=1; return
	end





%%%%%%%%%%%%%%%% SUBROUTINES

	function correctionStep()
			% Sub routine that calculates the correction step by
			% applying the pseudo-inverse matrix to the current
			% orbit deviation from the target-orbit 'R0'. The
			% correction step is then applied and the old RMS-beam
			% reading is prepended to 'BPMhist'.
			R = [B(1,:)'; B(2,:)'];
			R(isnan(R))=0;
			
			dphi = Mplus * ((R-par.R0).*par.weight);
		
			if par.scaleDisp~=0
				SC = SCsetCavs2SetPoints(SC,SC.ORD.Cavity,'Frequency',-par.scaleDisp*dphi(end),'add');
				dphi = dphi(1:end-1);
			end
			SC = SCsetCMs2SetPoints(SC,par.CMords{1},-dphi(1:length(par.CMords{1})      ),1,'add');
			SC = SCsetCMs2SetPoints(SC,par.CMords{2},-dphi(  length(par.CMords{1})+1:end),2,'add');
			BPMhist = circshift(BPMhist,1);
			BPMhist(1) = sqrt(mean(R.^2,1));
			ERROR = 0;
	end

	function res = isStable(n,eps)
		% Is true, if the coefficient of variation of the last 'n' RMS
		% BPM-readings is less that 'eps'.
		CV=var(BPMhist(1:n),1)/std(BPMhist(1:n));
		res = CV<eps;
	end

	function res = isConverged(n,eps)
		% Preliminary convergence criterion
		% TODO: This is not a very good convergence criterion. PhA
		CV=var(BPMhist(1:n),1)/std(BPMhist(1:n));
		res = CV<eps;
	end
end
