function SC = SCregisterSupport(SC,varargin)
% SCregisterSupport
% =================
%
% NAME
% ----
% SCregisterSupport - Register magnet support structures in SC
%
% SYNOPSIS
% --------
% `SC = SCregisterSupport(SC, type, ords [, options])`
%
%
% DESCRIPTION
% -----------
% Initializes magnet support structures such as sections, plinths and girders in SC. The function
% input be given as name-value pairs, starting with the structure type and structure ordinates
% defining start- end endpoints. Optional arguments are set as the uncertainties of e.g. girder
% offsets in the sigma structure `SC.SIG.Support`.
%
% INPUT
% -----
% `SC`::
% 	The SC base structure.
% `type`::
%   String specifying the support structure type. Valid are 'Plinth', 'Girder' or 'Section'.
% `ords`::
%   `[2xN]` array of ordinates defining start and end locations of `N` registered support structures.
%
%
% RETURN VALUE
% ------------
% `SC`::
% 	The base structure containing required information of support structure.
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'Offset'`::
%   A [1x3] array defining horizontal, vertical and longitudinal offset uncertainties for the start 
%   points or [2x3] array defining horizontal, vertical and longitudinal offset uncertainties for
%   the start end endpoints. If end points have dedicated uncertainties, *SCapplyErrors* applies 
%   random offset errors of both start end endpoints of the corresponding support structure, 
%   effectively  tilting the support structure. 
%   If only start points have asigned uncertainties, *SCapplyErrors* applies to the support
%   structure endpoints the same offset error as to the start points, resulting in a paraxial 
%   translation of the element. Only in this case dedicated `'Roll'` uncertainties may be given which
%   then tilt the structure around it's center.
%   The actual magnet or BPM offsets resulting from the support structure offsets is calculated in
%   *SCupdateSupport* by interpolating on a straight line between girder start- and endpoints. Note 
%   that the coordinate system change due to bending magnets are ignored in this calculation. Thus,
%   the accuracy of the result is limited if dipole magnets are involved. This may be particularily
%   true in case of large sections and/or longitudinal offsets.
% `'Roll'`::
%   [1x3] array [az,ax,ay] defineing roll (around z-axis), pitch (roll around x-axis) and yaw (roll
%   around y-axis) angle uncertainties.
%
% EXAMPLES
% --------
% Registers the girder start end endpoints defined in `ords` and assignes the horiatonal,
% vertical and longitudinal girder offset uncertainties `dX`, `dY` and `dZ`, respectively, to the 
% girder start points. When the support errors are applied the girder endpoints will get the same 
% offset error as the start points, resulting in a paraxial translation of the girder.
% ------------------------------------------------------------------
% SC = SCregisterSupport(SC,'Girder',ords,'Offset',[dX dY dZ]);
% ------------------------------------------------------------------
% Registers the section start- end endpoints defined in `ords` and assignes the horiatonal and
% vertical section offset uncertainties `dX` and `dY`, respectively, to the start points. When
% the support errors are applied the section endpoints will get the same offset as the start points.
% ------------------------------------------------------------------
% SC = SCregisterSupport(SC,'Section',ords,'Offset',[dX dY 0]);
% ------------------------------------------------------------------
% Registers the girder start end endpoints defined in `ords`, assignes the roll uncertainty `dPhi`
% and the horiatonal and vertical girder offset uncertainties `dX1` and `dY1`, respectively to the
% start points and `dX2` and `dY2` to the endpoints. When the support errors are applied, all
% girder start- and endpoints will get random offset errors and the resulting yaw and pitch angles 
% are calculated accordingly.
% ------------------------------------------------------------------
% SC = SCregisterSupport(SC,'Girder',ords,'Offset',[dX1 dY1 0; dX2 dY2 0],'Roll',[dPhi 0 0]);
% ------------------------------------------------------------------
% Registers the girder start end endpoints defined in `ords` and assignes the horiatonal,
% vertical and longitudinal girder offset uncertainties `dX`, `dY` and `dZ`, respectively, and the 
% roll, pitch and yaw angle uncertainties `az`, `ax` and `ay`. When the support errors are applied
% the girders will expirience a paraxial translation according to the offsets plus the proper
% rotations around the three x-, y- and z-axes.
% ------------------------------------------------------------------
% SC = SCregisterSupport(SC,'Girder',ords,'Offset',[dX dY dZ],'Roll',[az ax ay]);
% ------------------------------------------------------------------
%
%
% SEE ALSO
% --------
% *SCgetOrds*, *SCupdateSupport*, *SCgetSupportOffset*, *SCplotSupport*, *SCapplyErrors*, *SCregisterMagnets*, *SCgetTransformation*


	% Check for any input
	if length(varargin)<2
		return
	end

	checkInput()

	% Make sure that ordinates are within ring
	Nele = length(SC.RING);
	ords = mod(varargin{2}-1,Nele)+1;
	type = varargin{1};

	% Define start end end ordinates in SC
	SC.ORD.(type) = ords;

	% Loop over elements
	for ordPair=ords
		% Loop over start end end points
		for n=1:2
			% Initialize offset and roll field in lattice elements
			SC.RING{ordPair(n)}.([type 'Offset']) = [0 0 0]; % [x,y,z]
			SC.RING{ordPair(n)}.([type 'Roll'  ]) = [0 0 0]; % [az,ax,ay]
		end

		% Loop over input name/pair-values if given
		for i=3:2:(length(varargin)-1)
			% Check if outdated offset is provided
			if strcmp(varargin{i},'Offset') && size(varargin{i+1},2)==2
				warning('New error model requires ''Offset'' to be an array of [dx,dy,dz]. dz=0 added. This warning will be removed in future updates. Please update your code. Thanks!')
				varargin{i+1}(:,3) = 0; % TODO: double check!
			end
			% Define uncertainties for start points
			SC.SIG.Support{ordPair(1)}.([type varargin{i}]) = varargin{i+1}(1,:);
			% Check if endpoint uncertainties are given
			if size(varargin{i+1},1)==2
				% Define uncertainties for endpoints
				SC.SIG.Support{ordPair(2)}.([type varargin{i}]) = varargin{i+1}(2,:);
			end
		end
	end


	function checkInput()
		if ~any(strcmp(varargin{1},{'Girder','Plinth','Section'}))
			error('Unsupported structure type. Allowed are ''Girder'', ''Plinth'' and ''Section''.')
		end
		if isempty(varargin{2}) || size(varargin{2},1)~=2 
			error('Ordinates must be a 2xn array of ordinates.')
		end
		if mod(length(varargin),2)
			error('Optional input must be given as name-value pairs.')
		end
		
		if any(diff(varargin{2},1)<0)
			fprintf('%d ''%s'' endpoint(s) might be upstream of startpoint(s).',sum(diff(varargin{2},1)<0),varargin{1})
		end
		
		if any(strcmp(varargin,{'Offset'}))
			offset = varargin{find(strcmp(varargin,{'Offset'}))+1};
			if size(offset,2)~=3 || (size(offset,1)~=1 && size(offset,1)~=2)
				error('Support structure offset uncertainty of ''%s'' must be given as [1x3] (start end endpoints get same offset errors) or [2x3] (start end endpoints get independent offset errors) array.',varargin{1})
			end
		end
		
		if any(strcmp(varargin,{'Roll'}))
			if length(varargin{find(strcmp(varargin,{'Roll'}))+1})~=3
				error('''%s'' roll uncertainty must be [1x3] array [az,ax,ay] of roll (around z-axis), pitch (roll around x-axis) and yaw (roll around y-axis) angle.',varargin{1})
			end
		end
	end
end
