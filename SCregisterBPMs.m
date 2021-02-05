function SC = SCregisterBPMs(SC,BPMords,varargin)
% SCregisterBPMs
% ==============
%
% NAME
% ----
% SCregisterBPMs - Registers BPMs in SC
%
% SYNOPSIS
% --------
% `SC = SCregisterBPMs(SC,BPMords [, sigmas])`
%
%
% DESCRIPTION
% -----------
% Registers BPMs specified by `BPMords` in the `SC` structure and initializes all required fields 
% in the lattice elements. The ordinates of all registered BPMs are stored in `SC.ORD.BPM`. The 
% BPM realated fields in the lattice elements are
%
% `Noise`::
%   [1 x 2] array of hor./ver. turn-by-turn BPM noise
% `NoiseCO`::
%   [1 x 2] array of hor./ver. orbit BPM noise
% `CalError`::
%   [1 x 2] array of hor./ver. BPM calibration errors
% `Offset`::
%   [1 x 2] array of individual hor./ver. BPM offsets
% `SupportOffset`::
%   [1 x 2] array of hor./ver. BPM offsets which result from the corresponding girder offset at
%   the location of the BPMs, see *SCupdateSupport*.
% `Roll`::
%   BPM roll around z-axis w.r.t. the support structure
% `SupportRoll`::
%   BPM roll around z-axis which results from the corresponding support structure roll at
%   the location of the BPMs, see *SCupdateSupport*.
%
%
% INPUTS
% ------
% `SC`::       SC base structure.
% `BPMords`::  BPM ordinates in the lattice structure.
%
%
% UNCERTAINTIES
% -------------
% Additional name/vale-pairs are interpreted as uncertainties and passed to the sigma structure
% `SC.SIG` for the corresponding BPM. The function *SCapplyErrors* uses the fields of `SC.SIG` to
% randomly generate errors and applies them to the corresponding fields in `SC.RING`.
%
%
% RETURN VALUE
% ------------
% `SC`::
% 	The base structure containing required information of all BPMs.
%
%
% EXAMPLES
% --------
% Identify the ordinates of all elements named `BPM` and registers them as BPMs in `SC`
% ------------------------------------------------------------------
% ords = SCgetOrds(SC.RING,'BPM');
% SC = SCregisterBPMs(SC,ords);
% ------------------------------------------------------------------
%
% Register the BPMs specified in `ords` in `SC` and set the uncertanty of the offset to `500um` in 
% both planes. A subsequent call of *SCapplyErrors* would generate a random BPM offset errors with 
% `sigma=500um`.
% ------------------------------------------------------------------
% SC = SCregisterBPMs(SC,ords,'Offset',500E-6*[1 1]);
% ------------------------------------------------------------------
%   
% SEE ALSO
% --------
% *SCgetBPMreading*, *SCgetOrds*, *SCsanityCheck*, *SCapplyErrors*, *SCregisterSupport*, *SCupdateSupport*


	% Store BPM ordinates
	if isfield(SC,'ORD') && isfield(SC.ORD,'BPM')
		SC.ORD.BPM = sort(unique([SC.ORD.BPM BPMords]));
	else
		SC.ORD.BPM = BPMords;
	end
	
	% Loop over newly registered BPMs
	for ord = BPMords(:)'

		% Set name/pair-values in sigma structure
		if ~isempty(varargin)
			for i=1:2:(length(varargin)-1)
				SC.SIG.BPM{ord}.(varargin{i}) = varargin{i+1}(:)';
			end
		end
	
		% Set TBT BPM noise
		SC.RING{ord}.Noise = zeros(1,2);

		% Set orbit BPM noise
		SC.RING{ord}.NoiseCO = zeros(1,2);

		% Set BPM offset w.r.t. girder
		SC.RING{ord}.Offset = zeros(1,2);

		% Set girder offset at BPM location
		SC.RING{ord}.SupportOffset = zeros(1,2);

		% Set BPM roll w.r.t. girder
		SC.RING{ord}.Roll = 0;

		% Set girder roll at BPM location
		SC.RING{ord}.SupportRoll = 0;

		% Set BPM calibration
		SC.RING{ord}.CalError = zeros(1,2);
	end
	
end
