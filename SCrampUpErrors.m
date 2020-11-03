function [SC,ERROR] = SCrampUpErrors(SC,varargin)
% SCrampUpErrors
% ==============
%
% NAME
% ----
% SCrampUpErrors - Ramps up errors while applying trajectory/orbit feedback
%
% SYNOPSIS
% --------
% `[SC, ERROR] = SCrampUpErrors(SC [, options])`
%
%
% DESCRIPTION
% -----------
% This function may be used as a shortcut to the result of start-to-end correction chain. The
% input is the SC structure containing a `SC.RING` with already applied errors. *SCrampUpErrors*
% then ramps up these errors in steps while applying *SCfeedbackRun*. If ramping
% up fails, the number of steps is doubled and the procedure starts over until the number of steps
% exceeds 100.
%
%
% INPUT
% -----
% `SC`::
%	The base structure containing `SC.RING` with applied errors
%
% RETURN VALUES
% -------------
% `SC`::
%	The base structure containing `SC.RING` with applied errors and a corrected orbit
% `ERROR`::
%	Error value.
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'nStepsRamp'` (10)::
%	Number of steps for ramping up errors.
% `'alpha'` (10)::
%	Regularization parameter to get the pseudo inverse of the response matrix.
% `'target'` (0)::
%	(for feedback) break, if the RMS BPM reading reaches this value
% `'maxsteps'` (30)::
%	(for feedback) break, if this number of correction steps have been performed
% `'eps'` (`1e-5`)::
%	(for feedback) break, if the coefficient of variation of the RMS BPM reading is below
%	this value
% `'verbose'` (0)::
%	If true, debug information is printed.
%
%
% EXAMPLES
% --------
% Initialize `SC` with the lattice `RING`, apply the user defined function `registerRING` to
% register all required elements including uncertainties of the magnets, BPMs, RF, etc. In the next
% step the errors are applied by *SCapplyErrors*. Tracking is set to orbit mode and the previously
% applied errors are ramped up in 15 steps. If `err=0`, the output contains `SC.RING` with applied
% errors and a corrected closed orbit.
% ------------------------------------------------------------------
% SC = SCinit(RING);
% SC = registerRING(SC);
% SC = SCapplyErrors(SC);
% SC.INJ.trackMode = 'ORB';
% [SC,err] = SCrampUpErrors(SC,'nStepsRamp',15);
% ------------------------------------------------------------------
%
% SEE ALSO
% --------
% *SCregisterMagnets*, *SCregisterBPMs*, *SCregisterCAVs*, *SCapplyErrors*, *SCgetModelRM*, *SCgetPinv*,  *SCfeedbackRun*


	% Parse optional arguments
	p = inputParser;
	addOptional(p,'nStepsRamp',10);
	addOptional(p,'eps',1e-5);
	addOptional(p,'target',0);
	addOptional(p,'alpha',10);
	addOptional(p,'maxsteps',30);
	addOptional(p,'verbose',0);
	parse(p,varargin{:});
	par = p.Results;

	% Define fields which are to be ramped up
	errFieldsMag = {'CalErrorB','CalErrorA','RollAngle','T1','T2','PolynomAOffset','PolynomBOffset','MagnetOffset','MagnetRoll','SupportRoll','SupportOffset'};
	errFieldsBPM = {'Noise','NoiseCO','Offset','GirderOffset','Roll','CalError'};
	errFieldsRF  = {'Frequency','TimeLag','Voltage'};

	% Store SC structure with initial errors
	SC0 = SC;

	% Get model response matrix
	M = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'nTurns',SC.INJ.nTurns,'trackMode',SC.INJ.trackMode);

	% Calculate pseudo inverse
	Mplus = SCgetPinv(M,'alpha',par.alpha);

	% Loop over scaling steps
	for scale = linspace(1/par.nStepsRamp,1,par.nStepsRamp)

		% Print information
 		if par.verbose;fprintf('Ramping up errors with scaling factor %.2f.\n',scale);end


		% Scale magnet errors
		SC = scaleMagnets(SC,SC0,errFieldsMag,scale);

		% Scale BPM errors
		SC = scaleBPMs(SC,SC0,errFieldsBPM,scale);

		% Scale RF errors
		SC = scaleRF(SC,SC0,errFieldsRF,scale);

		% Scale injection errors
		SC = scaleInjection(SC,SC0,scale);

		% Scale circumference error
		SC = scaleCircumference(SC,SC0,scale);

		% Run feedback
		[CUR,ERROR] = SCfeedbackRun(SC,Mplus,'target',par.target,'maxsteps',par.maxsteps,'eps',par.eps,'verbose',par.verbose);
		if ~ERROR
			% Apply feedback result
			SC = CUR;
		else
			% Check if number of ramp-up steps exceeds reasonable values (due to recursive calling)
			if 2*par.nStepsRamp>100
				error('Ramping up failed at scaling %.2f with %d ramping steps. Try different feedback parameters.',scale,par.nStepsRamp)
			else
				fprintf('Feedback did not succeed at scaling %.2f. Trying with %d ramping steps.',scale,2*par.nStepsRamp)
				SC = SCrampUpErrors(SC0,varargin{:},'nStepsRamp',2*par.nStepsRamp);
				return
			end
		end
	end
end

% Scale magnet errors
function SC = scaleMagnets(SC,SC0,fields,scale)
	% Loop over magnets
	for ord=SC.ORD.Magnet
		% Loop over fields
		for field=fields
			if isfield(SC.RING{ord},field{1})
				SC.RING{ord}.(field{1}) = scale * SC0.RING{ord}.(field{1});
			end
		end
	end
	% Update magnets
	SC = SCupdateMagnets(SC);
end

% Scale BPM errors
function SC = scaleBPMs(SC,SC0,fields,scale)
	for field=fields
		SC.BPM.(field{1}) = scale * SC0.BPM.(field{1});
	end
end

% Scale RF errors
function SC = scaleRF(SC,SC0,fields,scale)
	for field=fields
		for ord=SC.ORD.Cavity
			SC.RING{ord}.(field{1}) = SC.IDEALRING{ord}.(field{1}) + scale * ( SC0.RING{ord}.(field{1}) - SC.IDEALRING{ord}.(field{1}) );
		end
	end
end

% Scale injection errors
function SC = scaleInjection(SC,SC0,scale)
	% Increase injectied beam initial trajectory
	SC.INJ.Z0 = SC0.INJ.Z0ideal + scale * (SC0.INJ.Z0 - SC0.INJ.Z0ideal);
	% Increase injectied beam jitter
	SC.INJ.randomInjectionZ = scale * SC0.INJ.randomInjectionZ;
	% Increase injectied beam size
	SC.INJ.beamSize = scale * SC0.INJ.beamSize;
end

% Scale circumference error
function SC = scaleCircumference(SC,SC0,scale)
	% Get circumferences
	D=0;D0=0;
	for ord=1:length(SC0.RING)
		D  = D  + SC0.RING{ord}.Length;
		D0 = D0 + SC0.IDEALRING{ord}.Length;
	end
	% Apply circumference error
	SC.RING = SCscaleCircumference(SC.RING,scale*(D-D0)+D0,'abs');
end
