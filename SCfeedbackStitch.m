function [SC,ERROR] = SCfeedbackStitch(SC,Mplus,varargin)
% SCfeedbackStitch
% ================
%
% NAME
% ----
% SCfeedbackStitch - achieves 2-turn transmission
%
% SYNOPSIS
% --------
% `[SC, ERROR] = SCfeedbackStitch(SC, Mplus [, options])`
%
%
% DESCRIPTION
% -----------
% The purpose of this function is to go from 1-turn transmission to 2-turn
% transmission. This is done by applying orbit correction based on the pseudo
% inverse trajectory response matrix 'Mplus' applied to the first BPMs in
% the 'SC.RING'. The reading of the BPMs in the second turn is corrected
% towards the reading of these BPMs in the first turn. This approach has been
% seen to be more stable than the direct application of the two-turn inverse
% response matrix to the two-turn BPM data.
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
% `'nBPMs'` (`4`)::
%	number of BPMs which are used in the second turn
% `'maxsteps'` (`30`)::
%	break, if this number of correction steps have been performed
% `'R0'` (zeros)::
%	target orbit in the format `[x_1 ... x_n y_1 ...y_n]`, where
%	`[x_i,y_i]` is the target position at the i-th BPM.
% `'nRepro'` (3)::
%	Number of consecutive beam injections for which the target should be reached
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
%	`maxsteps` was reached, without producing two-turn transmission.
%
% EXAMPLES
% --------
%
% Calculate the 2-turn response matrix and get the pseudo inverse using a Tikhonov regularization
% parameter of 10. Switch the injection pattern to 2 turns and apply the stitching using the first
% three BPMs, a maximum of 20 steps and print debug information.
% ------------------------------------------------------------------	
% RM2 = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'nTurns',2);
% Minv2 = SCgetPinv(RM2,'alpha',10);
% SC.INJ.nTurns = 2;
% [CUR,ERROR] = SCfeedbackStitch( SC,Minv2,'nBPMs',3,'maxsteps',20,'verbose',1);
% ------------------------------------------------------------------
%
% SEE ALSO
% --------
%  *SCgetPinv*, *SCgetModelRM*, *SCfeedbackRun*, *SCfeedbackFirstTurn*, *SCfeedbackBalance*, *SCgetBPMreading*, *SCsetCMs2SetPoints*


	% Parse optional arguments
	p = inputParser;
	addOptional(p,'R0',zeros(2,1));
	addOptional(p,'nBPMs',4);
	addOptional(p,'maxsteps',30);
	addOptional(p,'CMords',SC.ORD.CM);
	addOptional(p,'BPMords',SC.ORD.BPM);
	addOptional(p,'verbose',0);
	addOptional(p,'nRepro',3);
	parse(p,varargin{:});
	par=p.Results;


	if par.verbose; fprintf('SCfeedbackStitch: Start\n'); end;

	% Initialize error value
	ERROR = 1;

	% Initialize history of last-reached BPMs
	BPMhist = -1*ones(1,100);

	% Wiggle (see SCfeedbackFirstTurn) if initially there is not enough
	% transmission.
	B = SCgetBPMreading(SC,'BPMords',par.BPMords);
	if ~isSignal(B,par.nBPMs)
		if par.verbose; fprintf('SCfeedbackStitch: Wiggling\n'); end;
		wiggle();
	end

	% Fail if there is still not enough transmission
	B = SCgetBPMreading(SC,'BPMords',par.BPMords);
	if ~isSignal(B,par.nBPMs)
		if par.verbose; fprintf('SCfeedbackStitch: FAIL Wiggling failed\n'); end;
		ERROR=2; return;
	end


	% Main loop. In each step, a correction step is applied and the
	% break-conditions are checked.
	cnt=1;
	for steps=1:par.maxsteps
		B = SCgetBPMreading(SC,'BPMords',par.BPMords); % Inject...
		correctionStep(); % call correction subroutine.

		if isSetback(BPMhist)
			% Fail, is we have less transmission than before.
			if par.verbose; fprintf('SCfeedbackStitch: FAIL Setback\n'); end;
			ERROR=3; return;
		end

		if isRepro(BPMhist,par.nRepro) && isTransmit(BPMhist)
			% If we had full transmission in the last 'nRepro' injections
			% we succeed.
			if par.verbose; fprintf('SCfeedbackStitch: Success\n'); end;
			ERROR=0; return; % Success
		end

		cnt = cnt+1;
	end

	% If we have reached 'maxsteps' without producing full transmission,
	% that is considered a failure.
	if par.verbose; fprintf('SCfeedbackStitch: FAIL Reached maxsteps\n'); end;
	ERROR=1; return;



%%%%%%%%%%%%%%%% SUBROUTINES

	function correctionStep()
		% Sub routine that calculates the correction step by
		% applying the pseudo-inverse matrix to the current
		% first-turn orbit deviation from the target-orbit 'R0'
		% and second-turn deviation from the first turn.
		% The correction step is then applied and the old
		% RMS-beam reading is prepended to 'BRMShist'. The last
		BPMhist = logLastBPM(BPMhist,B);
		lBPM = size(B,2);
		Bx1 = B(1,1:lBPM/2);
		By1 = B(2,1:lBPM/2);
		Bx2 = B(1,(lBPM/2+1):end);
		By2 = B(2,(lBPM/2+1):end);
		DELTABx=Bx2-Bx1;
		DELTABy=By2-By1;
		DELTABx((par.nBPMs+1):end)=0;
		DELTABy((par.nBPMs+1):end)=0;
		R = [Bx1-par.R0(1,:),DELTABx,By1-par.R0(2,:),DELTABy]';
		R(isnan(R))=0;
		dphi = Mplus * R;
		SC = SCsetCMs2SetPoints(SC,par.CMords{1},-dphi(1:length(par.CMords{1})      ),1,'add');
		SC = SCsetCMs2SetPoints(SC,par.CMords{2},-dphi(  length(par.CMords{1})+1:end),2,'add');
	end


	function wiggle()
		% Wiggles a number of the last reached CMs around, based on a
		% golden-ratio-based covering of the 2D-torus.
		% TODO: Implement wiggle options as in `FirstTurn`
		pts = [[0;0] goldenDonut(50E-6,200E-6,32)];
		dpts = diff(pts,1,2);
		for nWiggleCM=[1 2 3 4 5 6 7 8]
			if par.verbose; fprintf('SCfeedbackStitch: Number of magnets used for wiggling: %d. \n',nWiggleCM); end
			CMords = getLastCMords(B,nWiggleCM);
			for i=1:size(dpts,2)
				for ord=CMords
					SC = SCsetCMs2SetPoints(SC,ord,dpts(1,i),1,'add');
					SC = SCsetCMs2SetPoints(SC,ord,dpts(2,i),2,'add');
				end
				W = SCgetBPMreading(SC,'BPMords',par.BPMords);
				BPMhist = logLastBPM(BPMhist,W);
				if isSignal(W,par.nBPMs) % TODO double check. Seems a bit iffy
					BPMhist = logLastBPM(BPMhist,SCgetBPMreading(SC,'BPMords',par.BPMords));
					BPMhist = logLastBPM(BPMhist,SCgetBPMreading(SC,'BPMords',par.BPMords));
					BPMhist = logLastBPM(BPMhist,SCgetBPMreading(SC,'BPMords',par.BPMords));
					if isRepro(BPMhist,3)
						BPMhist(1:3) = -1; % void last hist
						return
					end
				end
			end
		end
	end


	function ords = getLastCMords(B,n)
		% Returns last reached CMs
		dualCMords = intersect(par.CMords{:}); % Generate a list of CMs that can act in both planes.
		lastBPMidx = find(~isnan(B(1,:)),1,'last'); % Get index of last reached BPM
		if isempty(lastBPMidx) || lastBPMidx > length(par.BPMords)% If there is no beam loss in the first turn
			ords = dualCMords(end-n:end); ... just return the last n CMs
		else % in case the beam was lost in the first turn
			lastBPMord = par.BPMords(lastBPMidx); % We can then find the ordinate of the last BPM.
			lastCMidx  = find(dualCMords <= lastBPMord,1,'last'); % Find the last CM upstream of the last BPM.
			% Return the ``n'' upstream CMs (if there are so many)
			ords = dualCMords((lastCMidx-min(lastCMidx,n)+1):lastCMidx);
		end
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

function res = isNew(BPMhist)
	% True if, last correction step has reached a new BPM
	res = BPMhist(1)~=BPMhist(2);
end

function res = isSetback(BPMhist)
	% True if, last correction step resulted in less transmission
	res = BPMhist(1)~=0 && BPMhist(1)<BPMhist(3) && BPMhist(2)<BPMhist(3);
end

function res = isRepro(BPMhist,N)
	% True, if last-reached BPM is the same for N injections
	res = all(BPMhist(1:N)==BPMhist(1));
end

function res = isTransmit(BPMhist)
	% True, if we have full transmission
	res = BPMhist(1)==0;
end

function res = isSignal(B,nBPMs)
	% TODO
	lastBPMidx = find(~isnan(B(1,:)),1,'last');
	res = lastBPMidx >= size(B,2)/2 + nBPMs;
end

function BPMhist = logLastBPM(BPMhist,B)
	% Write last reached BPM to 'BPMhist'.
	% '0' indicates all BPMs have been reached.
	BPMhist = circshift(BPMhist,1);
	ord = getLastBPMord(B);
	if ord
		BPMhist(1)=ord;
	else
		BPMhist(1)=0;
	end
end

function ord = getLastBPMord(B)
	% Returns last-reached BPM
	ord = find(isnan(B(1,:)),1)-1;
end
