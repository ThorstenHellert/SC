function RING = setApertures_ALSU_SR(RING)
% setApertures_ALSU_SR
% ======================
%
% NAME
% ----
% setApertures_ALSU_SR - set apertures for the ALSU Storage Ring
%
% SYNOPSIS
% --------
% `RING = setApertures_ALSU_SR(RING)`
%
%
% DESCRIPTION
% -----------
% Defines the aperture (`EApertures`,`RApertures`) in the lattice structure for 
% the ALSU storage ring 
%
% RETURN VALUE
% ------------
% 'RING'::
% 	Contains lattice with aperture values.

	
	% Define physical aparture % % % % % % % % % % % % % % % % % % % % % %
	apArcVal = 6E-3;                 % Radius of aparture in arcs
	apIDVal  = 3E-3;                  % Radius of aparture in IDs
	apRecVal =  1E-3*[-2,6,-6,6];  % Rectangular aparture value of rectangluar aparture at thin septum [xmin,xmax,ymin,ymax]
	apCavVal = 6E-3;                 % Radius of aparture in Cavity straight section

	% Find all elements which have not 'IdentityPass'
	ordAll = setdiff(1:length(RING),findcells(RING,'PassMethod','IdentityPass'));

	% Find all straight sections
	ordStraights = sort(findcells(RING,'FamName','D11'));

	% Define septum magnet ordinates
	ordSeptum = ordStraights([1 2 end-1 end]);

	% Define ID ordinates
	ordIDs    = ordStraights(3:end-2);

	% Define cacvity straight section
	CAVords = findcells(RING,'Frequency');
	ordCavs = ordStraights([find(ordStraights<CAVords,4,'last'),find(ordStraights>CAVords,4,'first')]);

	% Define elements with elliptical aperture
	ordElliptical = setdiff(ordAll,[ordSeptum ordIDs]);

	% Set large elliptical apertures in lattice structure
	RING = setcellstruct(RING,'EApertures'  ,ordElliptical,	num2cell(repmat( [apArcVal apArcVal], length(ordElliptical),1) ,2)  ,1     ,1:2);

	% Set small elliptical apertures in lattice structure
	RING = setcellstruct(RING,'EApertures'  ,ordIDs,        num2cell(repmat( [apIDVal   apIDVal], length(ordIDs),1) ,2)  ,1     ,1:2);

	% Set elliptical apertures in cavity straight section in lattice structure
	RING = setcellstruct(RING,'EApertures'  ,ordCavs,       num2cell(repmat( [apCavVal apCavVal], length(ordCavs),1) ,2)  ,1     ,1:2);

	% Set rectangular apertures at septum in lattice structure
	RING = setcellstruct(RING,'RApertures'  ,ordSeptum,  	num2cell(repmat(apRecVal,length(ordSeptum),1) ,2)  ,1     ,1:4);

end
