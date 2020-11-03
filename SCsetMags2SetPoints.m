function SC = SCsetMags2SetPoints(SC,MAGords,type,order,setpoints,varargin)
% SCsetMags2SetPoints
% ===================
%
% NAME
% ----
% SCsetMags2SetPoints - Sets magnets to setpoints
%
% SYNOPSIS
% --------
% `SC = SCsetMags2SetPoints(SC, MAGords, type, order, setpoints [, options])`
%
%
% DESCRIPTION
% -----------
% Sets magnets (except CMs) as specified in `MAGords` to `setpoints` while `order` and `type` defines 
% which field entry should be used (see below). The setpoints may be given relative to their nominal 
% value or in absolute terms. If the considered quadrupole is a combined function magnet with
% non-zero bending angle and the kick compensation flag is switched on, the appropriate bending
% angle difference is calculated and the horizontal CM setpoint is changed accordingly to compensate
% for that dipole kick difference.
% If the setpoint of a skew quadrupole exceeds the limit specified in the corresponding lattice
% field `SkewQuadLimit`, the setpoint is clipped to that value and a warning is being printed (to
% switch off, use `warning('off','SC:SkewLim')`)
%
% INPUTS
% ------
%
% `SC`::         SC base structure.
% `MAGords`::    Magnet ordinates in the lattice structure.
% `order`::      Numeric value defining the order of the considered magnet: [1,2,3,...] => [dip,quad,sext,...]
% `type`::       Numeric value defining the type of the considered magnet: [1,2] => [skew/normal]
% `setpoints`::  Magnet setpoints, either single value applied to all magnets or [1 x N] array
%                matching the number of magnets specified in `MAGords`.
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'method'` (`'abs'`)::
%   Specifies how the setpoints should be applied. Possible are
%   * `'abs'`: Applies the setpoints specified in `setpoints`.
%   * `'rel'`: Applies the setpoints specified in `setpoints` relative to the nominal setpoints of
%          the considered magnets.
%   * `'add'`: Adds the setpoints specified in `setpoints` to the current setpoints of the
%          considered magnets.
% `'dipCompensation'` (0)::
%	Used for combined function magnets. If this flag is set and if there is a horizontal CM
%   registered in the considered magnet, the CM is used to compensate the bending angle difference
%   if the applied quadrupole setpoints differs from the design value.
%
%
% RETURN VALUES
% -------------
% `SC`::
%	The base structure containing lattice with modified and applied setpoints.
%
% EXAMPLES
% --------
%
% Identify the ordinates of all elements named `'SF'` and switch their sextupole component off.
% ------------------------------------------------------------------
% ords = SCgetOrds(SC.RING,'SF');
% SC = SCsetMags2SetPoints(SC,ords,2,3,0,'method','abs');
% ------------------------------------------------------------------
%
% Identify the ordinates of all elements named `QF` and `QD` and set their quadrupole component
% to 99% of their design value.
% ------------------------------------------------------------------
% ords = SCgetOrds(SC.RING,'QF|QD');
% SC = SCsetMags2SetPoints(SC,ords,2,2,0.99,'method','rel');
% ------------------------------------------------------------------
%
% SEE ALSO
% --------
% *SCupdateMagnets*, *SCregisterMagnets*, *SCsetCMs2SetPoints*


	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Input check

	% Parse input
	p = inputParser;
	addOptional(p,'method','abs');
	addOptional(p,'dipCompensation',0);
	parse(p,varargin{:});
	par = p.Results;

	inputCheck(SC,MAGords,type,order,setpoints,par)

	% Expand setpoints to all magnets
	if length(setpoints)==1
		setpoints = repmat(setpoints,1,length(MAGords));
	end

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Main function

	i=1;
	for ord=MAGords(:)'

		% Define nominal and current fields
		nomAB = [SC.RING{ord}.NomPolynomA',SC.RING{ord}.NomPolynomB'];
		curAB = [SC.RING{ord}.SetPointA'  ,SC.RING{ord}.SetPointB'  ];

		% Selct actual setpoint
		switch par.method
			case 'abs'
				setpoints(i) = setpoints(i);
			case 'rel'
				setpoints(i) = setpoints(i) * nomAB(order,type);
			case 'add'
				setpoints(i) = setpoints(i) + curAB(order,type);
			otherwise
				warning('Unsupported setpoint flag.\n')
		end

		% Check clipping of skew quads
		setpoints(i) = checkClipping(SC,ord,type,order,setpoints(i));

		% Compensate for bending kick difference.
		if par.dipCompensation
			SC = dipCompensation(SC,ord,setpoints(i));
		end

		% Apply setpoint.
		if type==1
			SC.RING{ord}.SetPointA(order) = setpoints(i);
		else
			SC.RING{ord}.SetPointB(order) = setpoints(i);
		end

		% Update magnets.
		SC = SCupdateMagnets(SC,ord);

		% Update loop index
		i = i + 1;
	end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions

% Input check
function inputCheck(SC,MAGords,type,order,setpoints,par)
	if order==1
		error('Function not specified for CMs. Use ''SCsetCMs2SetPoints'' instead.')
	end
	if par.dipCompensation && order~=2
		error('Dipole compensation only works for quadrupoles.')
	end
	for ord=MAGords(:)'
		if ~isfield(SC.RING{ord},'SetPointB')
			error('Lattice element %d is no registered magnet.',ord)
		end
		if length(SC.RING{ord}.SetPointB)<order
			error('Lattice element %d is no registered magnet with order %d.',ord,order)
		end
	end
	if length(setpoints)~=1 && length(setpoints)~=length(MAGords)
		error('''setpoints'' must be either a single value or match the length of magnet ordinates.')
	end
end

% Apply dipole compensation
function SC = dipCompensation(SC,ord,setpoint)
	% Check if diple compensation possible
	if ~(isfield(SC.RING{ord},'BendingAngle') && SC.RING{ord}.BendingAngle ~= 0 && ismember(ord,SC.ORD.CM{1}))
		return
	end

	% Calculate bending kick differnece for ideal magnet. See note-y18m08d20.
	idealKickDifference =  ( ( setpoint - ( SC.RING{ord}.SetPointB(2)-SC.RING{ord}.NomPolynomB(2) ) ) / SC.RING{ord}.NomPolynomB(2) - 1) * SC.RING{ord}.BendingAngle / SC.RING{ord}.Length;

	% Set dipole setpoint accordinly.
	[SC,~] = SCsetCMs2SetPoints(SC,ord, idealKickDifference*SC.RING{ord}.Length ,1,'add');
end

% Check for skew quadrupole clipping
function setpoint = checkClipping(SC,ord,type,order,setpoint)
	% If no skew quadrupole, return to main function
	if ~(type==1 && order==2)
		return
	end
	% Check if limit is specified and exceeded
	if isfield(SC.RING{ord},'SkewQuadLimit') && abs(setpoint)>abs(SC.RING{ord}.SkewQuadLimit)
		warning('SC:SkewLim','Skew quadrupole (ord: %d) is clipping',ord);
		setpoint = sign(setpoint) * SC.RING{ord}.SkewQuadLimit;
	end
end
