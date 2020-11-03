function SC = SCupdateCAVs(SC,varargin)
% SCupdateCAVs
% ============
%
% NAME
% ----
% SCupdateCAVs - Updates the cavities in the lattice
%
% SYNOPSIS
% --------
% `SC = SCupdateCAVs(SC [, ords])`
%
%
% DESCRIPTION
% -----------
% Updates the cavity fields `Voltage`, `Frequency` and `TimeLag` in `SC.RING` as specified in `ords`.
% If no ordinates are given explicitly, all registered cavities defined in `SC.ORD.Cavity` are 
% updated. For each cavity and each field, the setpoints, calibration errors and offsets are considered.
%
% INPUT
% -----
% `SC`::
% 	Base structure.
% `ords` (optional)::
% 	Cavity ordinates to be updated.
%
% RETURN VALUE
% ------------
% `SC`::
% 	Base structure with updated cavities.
%
% SEE ALSO
% --------
% *SCregisterCAVs*, *SCapplyErrors*


	% Check if ordinates are given explicitly. If not, take all registered cavities
	if isempty(varargin)
		ords = SC.ORD.Cavity;
	else
		ords = varargin{1};
	end

	% List of fields to be updated
	fields = {'Voltage','Frequency','TimeLag'};

	% Loop over all specified cavities
	for ord=ords
		for field=fields
			SC.RING{ord}.(field{1}) = SC.RING{ord}.([field{1} 'SetPoint']) * (1 + SC.RING{ord}.([field{1} 'CalError'])) + SC.RING{ord}.([field{1} 'Offset']);
		end
	end
end
