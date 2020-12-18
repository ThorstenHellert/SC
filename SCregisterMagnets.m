function SC = SCregisterMagnets(SC,MAGords,varargin)
% SCregisterMagnets
% =================
%
% NAME
% ----
% SCregisterMagnets - Register magnets in SC
%
% SYNOPSIS
% --------
% `SC = SCregisterMagnets(SC, MAGords [, options])`
%
%
% DESCRIPTION
% -----------
% Registers magnets specified by `MAGords` in the `SC` structure and initializes all required fields 
% in the lattice elements. The ordinates of all registered magnets are stored in `SC.ORD.Magnet`.
% The additional `SC` related fields in the lattice elements are
%
% `NomPolynomB`::
%   Nominal (design) `PolynomB` fields.
% `NomPolynomA`::
%   Nominal (design) `PolynomA` fields.
% `SetPointB`::
%   Setpoints for the `PolynomB` fields.
% `SetPointA`::
%   Setpoints for the `PolynomA` fields.
% `CalErrorB`::
%   Calibration error of the `PolynomB` fields wrt. the corresponding setpoints.
% `CalErrorA`::
%   Calibration error of the `PolynomA` fields wrt. the corresponding setpoints.
% `PolynomBOffset` (optional)::
%   Offset error of the `PolynomB` fields wrt. the corresponding setpoints.
% `PolynomAOffset` (optional)::
%   Offset error of the `PolynomA` fields wrt. the corresponding setpoints.
% `MagnetOffset`::
%   [1 x 2] array of horizontal and vertical magnet offsets (wrt. the support structure).
% `SupportOffset` (optional)::
%   [1 x 2] array of horizontal and vertical support structure offsets (if support structure is registered).
% `MagnetRoll`::
%   Magnet roll (wrt. the support structure).
% `SupportRoll` (optional)::
%   Support structure roll (if support structure is registered).
% `RollAngle`::
%   Magnet roll wrt. the design coordinate frame (used for field calculations).
% `CF` (optional)::
%   Flag identifying the corresponding magnet as a combined function dipole/quadrupole.
% `HCM` (optional)::
%   Flag identifying the corresponding magnet as a horizontal corrector magnet.
% `VCM` (optional)::
%   Flag identifying the corresponding magnet as a vertical corrector magnet.
% `SkewQuad` (optional)::
%   Flag identifying the corresponding magnet as a skew quadrupole corrector magnet.
% `MasterOf` (optional)::
%   Array of ordinates to which the corresponding magnet acts as master (split magnets).
% 
% Optional input arguments can be given as name/value-pairs and are either used
% to specify certain magnet types or to define uncertainties (see below).
% If CMs or skew quadrupole correctors are specified, the ordinates are also
% stored in the corresponding fields `SC.ORD.CM` and `SC.ORD.SkewQuad`,
% respectively.
%
%
% INPUTS
% ------
%
% `SC`::
%	SC base structure.
% `MAGords`::
%	Magnet ordinates in the lattice structure.
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'HCM'` ([])::
%	Magnet is identified as horizontal CM.
%	The corresponding value is the horizontal CM limit and stored in
%	`SC.RING{MAGords}.CMlimit(1)`. E.g. set limit to `Inf`.
% `'VCM'` ([])::
%	Magnet is identified as vertical CM.
%	The corresponding value is the vertical CM limit and stored in
%	`SC.RING{MAGords}.CMlimit(2)`. E.g. set limit to `Inf`.
% `'SkewQuad'` ([])::
%	Magnet is identified as skew quadrupole corrector.
%	The corresponding value is the skew quadrupole limit and stored in
%	`SC.RING{MAGords}.SkewLimit`. E.g. set limit to `Inf`.
% `'CF'` (0)::
%	If true, the magnet is identified as a combined function
%	dipole/quadrupole.
%	That implies that the bending angle depends on the quadrupole setpoint.
%	A variation from the design value will therefore result in a bending
%	angle error which is added to the `PolynomB(1)` field.
% `'MasterOf'` ([])::
%	The magnets `MAGords` are identified as a split magnets each with `N`
%	childs as specified in the corresponding value which must be a [`N` x
%	`length(MAGords)`] array.
%	The field calculation in *SCupdateMagnets* uses the setpoints and
%	errors of the master magnet to calculate the child fields.
%	The relative bending angle error of the master magnet e.g. is applied
%	on the corresponding child bending angle appropriately.
%	Split quadrupole magnets with different design gradients, however, can
%	currently not be updated correctly.
%
%
% UNCERTAINTIES
% -------------
% Every name/value-pair which is not explicitly mentioned in the options above
% is interpreted as an uncertainty and passed to the sigma structure `SC.SIG`
% for the corresponding magnets defined in `MAGords`.
% The function *SCapplyErrors* uses the fields of `SC.SIG` to randomly generate
% errors and applies them to the corresponding fields of the lattice elements.
%
%
% EXAMPLES
% --------
%
% Identify the ordinates of all elements named `QF` and register them in `SC`.
% ------------------------------------------------------------------
% ords = SCgetOrds(SC.RING,'QF');
% SC = SCregisterMagnets(SC,ords);
% ------------------------------------------------------------------
%
% Register the magnets specified in `ords` in `SC` and set the uncertainty of
% the quadrupole component to 1E-3.
% ----------------------------------------------------------------------------
% SC = SCregisterMagnets(SC,ords, ...
%	'CalErrorB',[0 1E-3]);
% ----------------------------------------------------------------------------
%
% Register split magnets.
% Identify the magnets named `BENDa` ([`1xN`] array `masterOrds`) and the
% magnets named `BENDb` and `BENDc` ([`2xN`] array `childOrds`) and register
% the `masterOrds` as the master magnets of the children in the corresponding
% columns of `childOrds`.
% The uncertanty of the bending angle is set to 1E-4.
% ----------------------------------------------------------------------------
% masterOrds = SCgetOrds(SC.RING,'BENDa');
% childOrds  = [SCgetOrds(SC.RING,'BENDb');SCgetOrds(SC.RING,'BENDc')];
% SC = SCregisterMagnets(SC,masterOrds, ...
%	'BendingAngle',1E-4,
%	'MasterOf',childOrds);
% ----------------------------------------------------------------------------
%
% Register the magnets specified in `ords` in `SC` and set the uncertainty of
% the quadrupole component to 1E-3 and the uncertainty of the bending angle to
% 1E-4.
% ----------------------------------------------------------------------------
% SC = SCregisterMagnets(SC,ords, ...
%	'CalErrorB',[0 1E-3], ...
%	'BendingAngle',1E-4);
% ----------------------------------------------------------------------------
%
% Register the magnets specified in `ords` in `SC` as combined function magnets
% and sets the uncertanty of the quadrupole component to 1E-3.
% ----------------------------------------------------------------------------
% SC = SCregisterMagnets(SC,ords, ...
%	'CF',1, ...
%	'CalErrorB',[0 1E-3]);
% ----------------------------------------------------------------------------
%
% Register the magnets specified in `ords` in `SC` and set the uncertanty of
% the skew quadrupole component to 2E-3 and the uncertanty of the sextupole
% component to 1E-3.
% ----------------------------------------------------------------------------
% SC = SCregisterMagnets(SC,ords, ...
%	'CalErrorA',[0 2E-3 0], ...
%	'CalErrorB',[0 0 1E-3]);
% ----------------------------------------------------------------------------
%
% Register the magnets specified in `ords` in `SC` as horizontal and vertical
% CMs, set their dipole uncertanties to 5% and 1%, respectively and define no
% CM limits.
% ----------------------------------------------------------------------------
% SC = SCregisterMagnets(SC,ords, ...
%	'HCM',Inf, ...
%	'VCM',Inf, ...
%	'CalErrorB',5E-2, ...
%	'CalErrorA',1E-2);
% ----------------------------------------------------------------------------
%
% Register the magnets specified in `ords` in `SC` as horizontal and vertical
% CMs, set their uncertanties to 5% and 1%, respectively and their limits to 1
% mrad. Furthermore, set the uncertanty of the skew quadrupole component to
% 2E-3 and the uncertanty of the sextupole component to 1E-3.
% ----------------------------------------------------------------------------
% SC = SCregisterMagnets(SC,ords, ...
%	'HCM',1E-3, ...
%	'VCM',1E-3, ...
%	'CalErrorB',[5E-2 0 1E-3], ...
%	'CalErrorA',[1E-2 2E-3 0]);
% ----------------------------------------------------------------------------
%
%
% SEE ALSO
% --------
% *SCgetOrds*, *SCupdateMagnets*, *SCsanityCheck*, *SCapplyErrors*,
% *SCregisterSupport*


	% Specify which optional arguments should not be written as uncertanties in SC.SIG
	keywords = {'HCM','VCM','CF','SkewQuad','MasterOf'};

	% Get name/value-pairs for sigma structure
	[nvpairs] = getSigmaPairs(keywords,varargin{:});

	% Loop over magnets
	for ord = MAGords(:)'

		% Set nominal polynom
		SC.RING{ord}.NomPolynomB = SC.RING{ord}.PolynomB;
		SC.RING{ord}.NomPolynomA = SC.RING{ord}.PolynomA;

		% Set setpoint
		SC.RING{ord}.SetPointB = SC.RING{ord}.PolynomB;
		SC.RING{ord}.SetPointA = SC.RING{ord}.PolynomA;

		% Set calibration factors
		SC.RING{ord}.CalErrorB = zeros(size(SC.RING{ord}.PolynomB));
		SC.RING{ord}.CalErrorA = zeros(size(SC.RING{ord}.PolynomA));

		% Set magnet and support offset
		SC.RING{ord}.MagnetOffset  = [0 0];
		SC.RING{ord}.SupportOffset = [0 0];

		% Set magnet, support and overall roll angle
		SC.RING{ord}.MagnetRoll  = 0;
		SC.RING{ord}.SupportRoll = 0;
		SC.RING{ord}.RollAngle   = 0;

		% Initialize T1 and T2 fields (not done in AT2.0 by default anymore)
		SC.RING{ord}.T1 = zeros(6,1);
		SC.RING{ord}.T2 = zeros(6,1);

		% Set optional attributes
		SC = setOptional(SC,ord,MAGords,varargin{:});

		% Set name/pair-values in sigma structure
		for i=1:2:(length(nvpairs)-1)
			SC.SIG.Mag{ord}.(nvpairs{i}) = nvpairs{i+1};
		end
	end

	% Store magnet ordinates
	SC = storeOrds(SC,MAGords(:)',varargin{:});
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions

% Sort options and sigma name/value-pair
function [nvpairs] = getSigmaPairs(keywords,varargin)
	nvpairs={};
	for n=1:2:length(varargin)
		if ~any(strcmp(varargin{n},keywords))
			% Write input argument in sigma name/value-pair
			nvpairs = horzcat(nvpairs,varargin{n},varargin{n+1});
		end
	end
end

% Set optional attributes to magnets
function SC = setOptional(SC,ord,MAGords,varargin)
	% Define combined function dipole
	if any(strcmp(varargin,'CF'))
		SC.RING{ord}.CombinedFunction = 1;
	end
	% Set CM limits if CMs are given
	if any(strcmp(varargin,'HCM'))
		SC.RING{ord}.CMlimit(1) = varargin{find(strcmp(varargin,'HCM'))+1};
	end
	if any(strcmp(varargin,'VCM'))
		SC.RING{ord}.CMlimit(2) = varargin{find(strcmp(varargin,'VCM'))+1};
	end
	% Set skew quad limits
	if any(strcmp(varargin,'SkewQuad'))
		SC.RING{ord}.SkewQuadLimit = varargin{find(strcmp(varargin,'SkewQuad'))+1};
	end
	% Check for master/child magnets
	if any(strcmp(varargin,'MasterOf'))
		% Fill master field in lattice structure with child ordinates
		SC.RING{ord}.MasterOf = varargin{find(strcmp(varargin,'MasterOf'))+1}(:,ord==MAGords)';
	end
end

% Store ordinates in SC structure
function SC = storeOrds(SC,MAGords,varargin)
	% Store magnet ordinates
	if isfield(SC,'ORD') && isfield(SC.ORD,'Magnet')
		SC.ORD.Magnet = sort(unique([SC.ORD.Magnet MAGords]));
	else
		SC.ORD.Magnet = MAGords;
	end

	% Store skew quadrupole ordinates
	if any(strcmp(varargin,'SkewQuad'))
		if isfield(SC,'ORD') && isfield(SC.ORD,'SkewQuad')
			SC.ORD.SkewQuad = sort(unique([SC.ORD.SkewQuad MAGords]));
		else
			SC.ORD.SkewQuad = MAGords;
		end
	end
	% Store horizontal CM ordinates
	if any(strcmp(varargin,'HCM'))
		if isfield(SC,'ORD') && isfield(SC.ORD,'CM')
			SC.ORD.CM{1} = sort(unique([SC.ORD.CM{1} MAGords]));
		else
			SC.ORD.CM{1} = MAGords;
		end
	end
	% Define vertical CM ordinates
	if any(strcmp(varargin,'VCM'))
		if isfield(SC,'ORD') && isfield(SC.ORD,'CM') && length(SC.ORD.CM)==2
			SC.ORD.CM{2} = sort(unique([SC.ORD.CM{2} MAGords]));
		else
			SC.ORD.CM{2} = MAGords;
		end
	end
end
