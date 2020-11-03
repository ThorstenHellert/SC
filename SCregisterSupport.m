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
%   A [2x1] array defining horizontal and vertical offset uncertainty for the start points or [2x2]
%   array defining horizontal and vertical offset uncertainties for the start end endpoints. If end
%   points have dedicated uncertainties, *SCapplyErrors* applies random offset errors of both start
%   end endpoints of the corresponding support structure, effectively giving the support structure a
%   pseudo-pitch. If only start points have asigned uncertainties, the support structure will get a
%   pure offset.
% `'Roll'`::
%   Single value defining roll uncertainty. Note that currently only girder rolls are implemented.
%
%
% EXAMPLES
% --------
% Registers the girder start end endpoints defined in `ords`.
% ------------------------------------------------------------------
% SC = SCregisterSupport(SC,'Girder',ords);
% ------------------------------------------------------------------
% Registers the section start- end endpoints defined in `ords` and assignes the horiatonal and
% vertical section offset uncertainties `dX` and `dY`, respectively, to the start points. When
% the support errors gets applied the section endpoints will get the same offset as the start points.
% ------------------------------------------------------------------
% SC = SCregisterSupport(SC,'Section',ords,'Offset',[dX dY]);
% ------------------------------------------------------------------
% Registers the girder start end endpoints defined in `ords`, assignes the roll uncertainty `dPhi`
% and the horiatonal and vertical girder offset uncertainties `dX1` and `dY1`, respectively to the
% start points and `dX2` and `dY2` to the endpoints. When the support errors gets applied, all
% girder start- and endpoints will get random offset errors, hence a pseudo girder pitch is applied.
% ------------------------------------------------------------------
% SC = SCregisterSupport(SC,'Girder',ords,'Offset',[dX1 dY1 ; dX2 dY2],'Roll',dPhi);
% ------------------------------------------------------------------
%
%
% SEE ALSO
% --------
% *SCgetOrds*, *SCgetSupportOffset*, *SCplotSupport*, *SCapplyErrors*, *SCregisterMagnets*


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
			SC.RING{ordPair(n)}.([type 'Offset']) = [0 0];
			SC.RING{ordPair(n)}.([type 'Roll'  ]) = 0;
		end

		% Loop over input name/pair-values if given
		for i=3:2:(length(varargin)-1)
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
		if size(varargin{2},1)~=2
			error('Ordinals must be a 2xn array of ordinates.')
		end
		if mod(length(varargin),2)
			error('Optional input must be given as name-value pairs.')
		end
		if any(strcmp(varargin,{'Offset'}))
			offset = varargin{find(strcmp(varargin,{'Offset'}))+1};
			if size(offset,2)~=2 || (size(offset,1)~=1 && size(offset,1)~=2)
				error('Support structure offset uncertainty must be given as [1x2] or [2x2] array.')
			end
		end
		if any(strcmp(varargin,{'Roll'})) && length(varargin{find(strcmp(varargin,{'Roll'}))+1})~=1
			error('Support structure roll uncertainty must be single value.')
		end
	end
end
