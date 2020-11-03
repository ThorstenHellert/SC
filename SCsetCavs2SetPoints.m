function SC = SCsetCavs2SetPoints(SC,CAVords,type,setpoints,varargin)
% SCsetCavs2SetPoints
% ===================
%
% NAME
% ----
% SCsetCavs2SetPoints - Set RF properties to setpoints
%
% SYNOPSIS
% --------
% `SC = SCsetCavs2SetPoints(SC, CAVords, type, setpoints [, mode])`
%
%
% DESCRIPTION
% -----------
% Set the setpoints of `Voltage`, `Frequency` or `TimeLag` as specified in `'type'` of the rf
% cavities specified in `CAVords`. If only a single setpoint is given for multiple cavities,
% the setpoint is applied to all cavities.
%
% INPUTS
% ------
% `SC`::         SC base structure
% `CAVords`::    Array of cavity ordinates in the lattice structure
% `type`::       String ('Voltage', 'Frequency' or 'TimeLag') specifying which cavity field should be set.
% `setpoints`::  Setpoints (array or single value for all cavities)
%
% MODE
% ----
% `'abs'` (default)::
%   Use absolute setpoint
% `'rel'`::
%   Use relative setpoint to current value
% `'add'`::
%   Add setpoints to current value
%
% RETURN VALUE
% ------------
% `SC`::
% 	The modified SC structure.
%
%
% EXAMPLES
% --------
% Sets the time lag of all cavities registered in SC to zero.
% -----------------------
% SC = SCsetCavs2SetPoints(SC, SC.ORD.Cavity, 'TimeLag', 0);
% -----------------------
%
% Adds 1kHz to the frequency of the first cavity.
% -----------------------
% SC = SCsetCavs2SetPoints(SC, SC.ORD.Cavity(1), 'Frequency', 1E3, 'add');
% -----------------------

	% Use absolute setpoints if not specified differently
	if isempty(varargin)
		mode = 'abs';
	else
		mode = varargin{1};
	end

	% If only single setpoint is given, use setpoint for all ordinates
	if length(setpoints)==1
		setpoints = repmat(setpoints,1,length(CAVords));
	end

	% Set setpoints according to typoe and options
	i = 1;
	for ord = CAVords
		switch mode
			case 'abs'
				SC.RING{ord}.([type 'SetPoint']) = setpoints(i);
			case 'rel'
				SC.RING{ord}.([type 'SetPoint']) = setpoints(i) * SC.RING{ord}.([type 'SetPoint']);
			case 'add'
				SC.RING{ord}.([type 'SetPoint']) = setpoints(i) + SC.RING{ord}.([type 'SetPoint']);
			otherwise
				warning('Unsupported setpoint type.\n')
		end
		i = i + 1;
	end

	% Update cavity fields
	SC = SCupdateCAVs(SC,CAVords);
	
end
