function [SC,wOrds] = IDLib_closeID(SC,ID,varargin)
% IDLib_closeID
% =============
%
% NAME
% ----
% IDLib_closeID - 'Closes' the ID
%
% SYNOPSIS
% --------
% `[SC,wOrds] = IDLib_closeID(SC,ID [,options]))`
%
% DESCRIPTION
% -----------
% This function 'closes' the ID by either setting a kickmap or changing the passmethod for SBENDs.
%
% INPUT
% -----
% `SC`::
% 	The SC base structure.
% `ID`::
% 	String with ID name.
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `IDstrength` (1)::
%	Scaling factor for ID (only for series of SBENDs)
% `'IDmap'` (`''`)::
%	Name of kick map to be loaded.
% `'PassMethod'` (`'ThinEPU2Pass'`)::
%	Passmethod for kickmaps.
%
% RETURN VALUE
% ------------
% `SC`::
% 	The SC base structure with 'closed' ID.
%
% SEE ALSO
% --------
% *IDlib_includeIDs*	

	% Parse optional arguments
	p = inputParser;
	addOptional(p,'IDstrength',1);
	addOptional(p,'IDmap','');
	addOptional(p,'PassMethod','ThinEPU2Pass');
	parse(p,varargin{:});
	par = p.Results;
	
	
	if strncmp(ID,'LEDA',4) || strcmp(ID,'COSMIC') || strcmp(ID,'TENDER') || strcmp(ID,'EPU36') || strcmp(ID,'XType')
		% Import and convert data
		map  = importdata(par.IDmap);
		brho = 6.67;
		nx = map.data(6,1);
		ny = map.data(8,1);
		[XEPU, YEPU] = meshgrid(map.data(11,1+(1:nx)),map.data(11+(1:ny),1));
		PXEPUS = map.data(12:12+ny-1,1+1:1+nx)/brho/brho;
		PYEPUS = map.data(12+ny+3:12+ny+3+ny-1,1+1:1+nx)/brho/brho;
		
		% Create lattice element
		EPUKICK6 = epukick(ID, nx, ny, XEPU, YEPU, PXEPUS, PYEPUS, par.PassMethod);
		global FAMLIST
		EPU = FAMLIST{EPUKICK6}.ElemData;
		
		% Get ID ordinates
		wOrds = SCgetOrds(SC.RING,ID);
		
		% Insert ID
		for ord=wOrds
			SC.RING{ord} = EPU;
		end
		
	elseif strcmp(ID,'EPU35') || strcmp(ID,'EPU90') || strcmp(ID,'EPU50') || strcmp(ID,'EPU70')
		% Import and convert data
		map = importdata(par.IDmap);
		nx  = length(unique(map(:,1)));
		ny  = length(unique(map(:,2)));
		[XEPU, YEPU] = meshgrid(1E-3*sort(unique(map(:,1))),1E-3*sort(unique(map(:,2))));
		imode = 1;
		nmode = nx*ny;
		if YEPU(1,1)<YEPU(end,1)
			YEPU = -YEPU;
		end
		for n=1:size(map,1)
			nX = find(1E-3*map(n,1)==XEPU(1,:));
			nY = find(1E-3*map(n,2)==YEPU(:,1));
			if isempty(nX) || isempty(nY)
				warning('Something went wrong at index %d.',n)
			end
			PXEPUS(nY,nX) = map(n,6);
			PYEPUS(nY,nX) = map(n,8);
		end
		
		% Create lattice element
		EPUKICK6 = epukick(ID, nx, ny, XEPU, YEPU, PXEPUS, PYEPUS, par.PassMethod);
		global FAMLIST
		EPU = FAMLIST{EPUKICK6}.ElemData;
		
		% Get ID ordinates
		wOrds = SCgetOrds(SC.RING,ID);
		
		% Insert ID
		for ord=wOrds
			SC.RING{ord} = EPU;
		end
	else
		% Get ID ordinates
		wOrds = SCgetOrds(SC.RING,[ID '$']);
		
		% Close ID
		for ord=wOrds
			SC.RING{  ord}.Class      = 'Bend';
			SC.RING{  ord}.PassMethod = 'BndMPoleSymplectic4RadPass';
			
			SC.RING{ord}.BendingAngle  = par.IDstrength * SC.RING{ord}.BendingAngle;
			SC.RING{ord}.EntranceAngle = par.IDstrength * SC.RING{ord}.EntranceAngle;
			SC.RING{ord}.ExitAngle     = par.IDstrength * SC.RING{ord}.ExitAngle;
		end
	end
end
