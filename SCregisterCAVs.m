function SC = SCregisterCAVs(SC,CAVords,varargin)
% SCregisterCAVs
% ==============
%
% NAME
% ----
% SCregisterCAVs - Register cavities in SC
%
% SYNOPSIS
% --------
% `SC = SCregisterCAVs(SC,CAVords [, sigmas])`
%
%
% DESCRIPTION
% -----------
% Register cavities specified in `CAVords` in `SC` by initializing all required fields in the 
% corresponding cavity lattice elements and storing the ordinates in `SC.ORD.Cavity`. The additional 
% fields in the lattice elements are
%
% `VoltageSetPoint`::
%   Setpoint of cavity voltage
% `VoltageOffset`::
%   Offset of cavity voltage wrt. to the setpoint
% `VoltageCalError`::
%   Calibration error of cavity voltage wrt. to the setpoint
% `FrequencySetPoint`::
%   Setpoint of cavity frequency
% `FrequencyOffset`::
%   Offset of cavity frequency wrt. to the setpoint
% `FrequencyCalError`::
%   Calibration error of cavity frequency wrt. to the setpoint
% `TimeLagSetPoint`::
%   Setpoint of cavity phase (`TimeLag`)
% `TimeLagOffset`::
%   Offset of cavity phase wrt. to the setpoint
% `TimeLagCalError`::
%   Calibration error of cavity phase wrt. to the setpoint
%
%
% INPUTS
% ------
% `SC`::       SC base structure.
% `CAVords`::  Cavity ordinates in the lattice structure.
%
%
% UNCERTAINTIES
% -------------
% Additional name/vale-pairs are interpreted as uncertainties and passed to the sigma structure
% `SC.SIG` for the corresponding cavity. The function *SCapplyErrors* uses the fields of `SC.SIG`
% to randomly generate errors and applies them to the corresponding fields of the lattice element.
% By default a 2 sigma cutoff is applied. The user can specify a different cutoff by giving the 
% uncertainty as a cell structure, e.g. {deltaF,nSig}, with nSig being the cutoff (see examples 
% below).
%
%
% RETURN VALUE
% ------------
% `SC`::
% 	The base structure containing required information of all cavities.
%
%
% EXAMPLES
% --------
% Identify the ordinates of all elements named `'CAV'` and register them as cavities in `SC`
% ------------------------------------------------------------------
% ords = SCgetOrds(SC.RING,'CAV');
% SC = SCregisterCAVs(SC,ords);
% ------------------------------------------------------------------
%
% Register the cavities specified in `ords` in `SC` and sets the uncertanty of the frequency offset 
% to 1kHz. A subsequent call of *SCapplyErrors* would generate a random frequncy offset error with 
% `sigma=1kHz`.
% ------------------------------------------------------------------
% SC = SCregisterCAVs(SC,ords,'FrequencyOffset',1E3);
% ------------------------------------------------------------------
%
% Register the cavities specified in `ords` in `SC` and sets the uncertanty of the frequency offset 
% to 1kHz. A subsequent call of *SCapplyErrors* would generate a random frequncy offset error with 
% `sigma=1kHz` and a 3 sigma cutoff.
% ------------------------------------------------------------------
% SC = SCregisterCAVs(SC,ords,'FrequencyOffset',{1E3,3});
% ------------------------------------------------------------------
%
% Register the cavities specified in `ords` in `SC` and sets the uncertanty of the timelag offset 
% to 0.3m. A subsequent call of *SCapplyErrors* would generate a random timelag offset error 
% ('phase error') with `sigma=0.3m` and a 3 sigma cutoff.
% ------------------------------------------------------------------
% SC = SCregisterCAVs(SC,ords,'TimeLagOffset',{0.3,3});
% ------------------------------------------------------------------
%
% SEE ALSO
% --------
% *SCgetOrds*, *SCsanityCheck*, *SCapplyErrors*



	% Store cavity ordinates
	if isfield(SC,'ORD') && isfield(SC.ORD,'Cavity')
		SC.ORD.Cavity = sort(unique([SC.ORD.Cavity CAVords]));
	else
		SC.ORD.Cavity = sort(CAVords);
	end

	% List of fields for lattice element
	fields = {'Voltage','Frequency','TimeLag'};

	% Loop over newly registered cavities
	for ord = CAVords(:)'
		% Initialize required fields in lattice structure
		for field=fields
			SC.RING{ord}.([field{1} 'SetPoint']) = SC.RING{ord}.(field{1});
			SC.RING{ord}.([field{1} 'Offset'  ]) = 0;
			SC.RING{ord}.([field{1} 'CalError']) = 0;
		end

		% Set name/pair-values in sigma structure
		if ~isempty(varargin)
			for i=1:2:(length(varargin)-1)
				SC.SIG.RF{ord}.(varargin{i}) = varargin{i+1}(:)';
			end
		end
	end
	
end % End of function

