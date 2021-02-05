function RING = defineApertures_ALSU_AR(RING)
% Defines the aperture (`EApertures`,`RApertures`) in the lattice structure for
% the ALSU accumulator ring (AR)

	
	% Define physical aparture %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	apArcVal      = 12.6E-3*[1 1];     % Circular aparture in arcs
	apStraightVal = 12.6E-3*[1 1];     % Circular aparture in straights
	apDIPVal      = [18.3E-3 5.6E-3];  % Eliptical aparture in dipoles
	injAp         = 1E-3*[20.7, 5.1];  % Eliptical aparture in injection section
	
	apSeptumVal   =  1E-3*[-20,8,-15,15];   % Rectangular aparture value of rectangluar aparture at thin septum [xmin,xmax,ymin,ymax]  
	apKickerVal   =  1E-3*[-20,20,-3,3];    % Rectangular aparture value of rectangluar aparture at injection/extraction kickers [xmin,xmax,ymin,ymax]
	
	
	% Identify ordinates of elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Find all elements which have not 'IdentityPass'
	ordAll = setdiff(1:length(RING),findcells(RING,'PassMethod','IdentityPass'));
	
	% Define straight section ordinates
	ordStraights = SCgetOrds(RING,'L1');
	
	% Define inj./extr. kicker ordinates
	ordKicker   = SCgetOrds(RING,'KICKER|preDK|^DK|Vkick');
	
	% Define septum ordinates
	ordSeptum = SCgetOrds(RING,'SEPTUM');
	
	% Define dipole ordinates
	DIPords = SCgetOrds(RING,'BEND');
	
	% Define ordinates between injection septum and dipole
	tmpINJ = SCgetOrds(RING,'INJ');	
	injOrds = tmpINJ:(DIPords(find(DIPords>tmpINJ,1))-1);
	% Define ordinates between 1st and 2nd dipole
	tmp = DIPords(find(DIPords>tmpINJ,2));
	injOrds2 = (tmp(1)+1):(tmp(2)-1);
	
	% Define elements with elliptical aperture
	ordArc = setdiff(ordAll,[DIPords(:);ordKicker(:);ordStraights(:);ordSeptum(:)]);

	
	% Set apertures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Set straight section elliptical apertures in lattice structure
	for ord=ordStraights
		RING{ord}.EApertures = apStraightVal;
	end
	
	% Set arc elliptical apertures in lattice structure
	for ord=ordArc
		RING{ord}.EApertures = apArcVal;
	end
	
	% Set dipole elliptical apertures in lattice structure
	for ord=DIPords
		RING{ord}.EApertures = apDIPVal;
	end
	
	% Set dipole radiation shield in lattice structure
	for ord=(DIPords+1)
		RING{ord}.EApertures = apArcVal - [1.7E-3 0];
	end
	
	% Set special injection apertures between INJ and DK
	for ord=injOrds(1):injOrds(end)
		if isfield(RING{ord},'RApertures')
			RING{ord} = rmfield(RING{ord},'RApertures');
		end
		RING{ord}.EApertures = injAp;
	end
	
	% Set special injection apertures between first two dipoles
	for ord=injOrds2
		if isfield(RING{ord},'RApertures')
			RING{ord} = rmfield(RING{ord},'RApertures');
		end
		RING{ord}.EApertures = apArcVal + [2E-3 0];
	end
	
	% Set rectangular apertures at inj./extr. kicker in lattice structure
	for ord=ordKicker
		RING{ord}.RApertures = apKickerVal;
	end
	
	% Set rectangular apertures at septum in lattice structure
	for ord=ordSeptum
		RING{ord}.RApertures = apSeptumVal;
	end
	
end
