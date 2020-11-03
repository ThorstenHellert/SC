function [SC,setpoints] = SCsetCMs2SetPoints(SC,CMords,setpoints,nDim,varargin)
% SCsetCMs2SetPoints
% ==================
%
% NAME
% ----
% SCsetCMs2SetPoints - Sets dipole corrector magnets to different setpoints
%
% SYNOPSIS
% --------
% `SC = SCsetCMs2SetPoints(SC, CMords, setpoints, nDim [, mode])`
%
%
% DESCRIPTION
% -----------
% Sets horizontal or vertical CMs as specified in `CMords` and `nDim`, respectively, to `setpoints` 
% [rad] and updates the magnetic fields. If the corresponding setpoint exceeds the CM limit 
% specified in the corresponding lattice field `CMlimit`, the CM is clipped to that value
% and a warning is being printed (to switch off, use `warning('off','SC:CM1'))`. Positive setpoints
% will results in kicks in the positive horizontal or vertical direction.
%
%
% INPUTS
% ------
% `SC`::         SC base structure
% `CMords`::     Array of CM ordinates in the lattice structure
% `setpoints`::  CM setpoints (array or single value for all CMs) [rad]
% `nDim`::       Integer specifying CM dimension ([1|2] -> [hor|ver])
%
%
% RETURN VALUES
% -------------
% `SC`::
%	The lattice structure with modified and applied setpoints
% `setpoints`::
%	The list of acutal setpoints applied to the magnets after possible clipping [rad]
%
%
% MODE
% ----
% `'abs'` (default)::
%   Use absolute setpoints
% `'rel'`::
%   Use setpoints relative to current value
% `'add'`::
%   Add setpoints to current value
%
%
% EXAMPLES
% --------
% Set all registered horizontal CMs to zero.
% ------------------------------------------------------------------
% SC = SCsetCMs2SetPoints(SC,SC.ORD.CM{1},0,1);
% ------------------------------------------------------------------
%
% Add 10urad to the fourth registered vertical CM.
% ------------------------------------------------------------------
% SC = SCsetCMs2SetPoints(SC,SC.ORD.CM{2}(4),1E-5, 2,'add');
% ------------------------------------------------------------------
%
%
% SEE ALSO
% --------
% *SCregisterMagnets*, *SCsetMags2SetPoints*, *SCupdateMagnets*, *SCgetCMSetPoints*


	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Input check

	% If only single setpoint is given, use setpoint for all CMs
	if length(setpoints)==1
		setpoints = repmat(setpoints,size(CMords));
	end

	% Set setpoint flag
	if isempty(varargin)
		method = 'abs';
	else
		method = varargin{1};
	end

	i=1;
	for ord=CMords(:)'

		% Get current fields
		curAB = [SC.RING{ord}.SetPointB'  ,SC.RING{ord}.SetPointA'  ];

		% Get setpoint normalization factor (field <-> kick)
		if strcmp(SC.RING{ord}.PassMethod,'CorrectorPass')
			normBy = [1 1];
		else
			normBy = [-1 1]*SC.RING{ord}.Length; % positive setpoint -> positive kick -> negative horizontal field
		end

		% Select actual setpoint based on setpoint flag
		switch method
			case 'abs'
				setpoints(i) = setpoints(i);
			case 'rel'
				setpoints(i) = setpoints(i) * curAB(1,nDim) * normBy(nDim);
			case 'add'
				setpoints(i) = setpoints(i) + curAB(1,nDim) * normBy(nDim);
			otherwise
				error('Unsupported method: ''%s''. Allowed are ''abs'',''rel'' and ''add''.',method)
		end

		% Check clipping
		if isfield(SC.RING{ord},'CMlimit') && abs(setpoints(i))>abs(SC.RING{ord}.CMlimit(nDim))
			warning('SC:CM1','CM (ord: %d / dim: %d) is clipping',ord,nDim);
			setpoints(i) = sign(setpoints(i)) * SC.RING{ord}.CMlimit(nDim);
		end

		% Apply setpoint.
		if nDim==1
			SC.RING{ord}.SetPointB(1) = setpoints(i) / normBy(nDim);
		else
			SC.RING{ord}.SetPointA(1) = setpoints(i) / normBy(nDim);
		end

		% Update loop index
		i = i + 1;
	end

	% Update magnets
	SC = SCupdateMagnets(SC,CMords);
	
end
