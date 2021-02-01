function RING = SCsetMultipoles_ALSU_AR(RING,ords,AB,varargin)
%
% ==> The function *SCsetMultipoles* has been updated in 21/02 to allow for dynamically updated systematic 
% ==> multipole errors. In order to keep the ALSU-AR example consistent with the PRAB paper we decided
% ==> to stay with the old higer order mulitpole error model.
%
% SCsetMultipoles_ALSU_AR
% =======================
%
% NAME
% ----
% SCsetMultipoles_ALSU_AR - sets multipole errors in lattice elements as used in ALSU-AR example/paper
%
% SYNOPSIS
% --------
% `RING = SCsetMultipoles_ALSU_AR(RING, ords, AB [, options])`
%
%
% DESCRIPTION
% -----------
% Applies multipole errors specified in `AB` in the lattice elements `ords` of `RING` depending on
% the specified options.
%
%
% INPUTS
% ------
% `RING`::  Lattice cell structure.
% `ords`::  Ordinates of the considered magnets.
% `AB`::    [N x 2] array of multipole errors.
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'method'` (`'scaleRnd'`)::
%   Specifies which multipole error method should be applied. Possible are
%   * `'relToNom'`: Applies the multipoles relative to the design primary field component defined by option
%                   `'order'` and `'type'`. The multipole table entry of the primary component is set to zero.
%                   This option is intended for normalized systematic multipole error tables.
%   * `'sysRnd'`:   Scales the multipoles with a random number within the defined range specified by
% 			        `'scale'` and zeros the primary field component defined by option `'order'`
%                   and `'type'`. This option is intended for normalized systematic multipole error tables
%                   and accounts, e.g. for systematic multipole errors of CMs or skew quadrupoles without
%                   using their actual setpoint to calculate the multipoles. It is assumed that a rms
%                   value of the CM setpoints is sufficient to describe the effect of CM multipole errors
%                   on the lattice.
%   * `'scaleRnd'`: Randomly scales each entry of the multipole error table within the defined range
%                   specified by `'scale'`. This option is intended for random multipole errors of magnets.
% `'scale'` ([])::
%   Defines the scaling range of the multipoles.
% `'order'` ([])::
%   Numeric value defining the order of the considered magnet: [1,2,3,...] => [dip,quad,sext,...]
% `'type'` ([])::
%   Numeric value defining the type of the considered magnet: [1,2] => [normal/skew]
%
%
% SEE ALSO
% --------
% *SCmultipolesRead*, *SCupdateMagnets*


	p = inputParser;
	addParameter(p,'method','scaleRnd');
	addOptional(p,'scale',[]);
	addOptional(p,'order',[]);
	addOptional(p,'type',[]);
	parse(p,varargin{:});
	par = p.Results;

	checkInput(par);

	% Loop over magnets
	for ord=ords
		% Check which multipole error method should be applied
		switch p.Results.method
			% Relative to the specified main component
			case 'relToNom'

				% Identify coefficient for scaling HOMs
				if par.type==1 % Skew
					coeff = RING{ord}.NomPolynomA(par.order);
				elseif par.type==2 % Normal
					coeff = RING{ord}.NomPolynomB(par.order);
				end

				% Scale entries according to nominal coefficent
				addAB = coeff * AB;

				% Set nominal coefficent to zero (we dont want to add to it)
				addAB(par.order,par.type) = 0;

			case 'sysRnd'

				% Scale all HOMs with the same random number generated within user supplied range
				addAB = par.scale * SCrandnc(2,1,1) * AB;

				% Set nominal coefficent to zero (we dont want to add to it)
				addAB(par.order,par.type) = 0;

			case 'scaleRnd'

				% Scale all HOMs randomly within user supplied range
				addAB = par.scale * SCrandnc(2,size(AB)) .* AB;

		end
		% Apply multipole errros to lattice
		RING = applyMultipoles(RING,ord,addAB);
	end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions

% Check if input looks reasonable
function checkInput(par)
	switch par.method
		case 'relToNom'
			if isempty(par.order) || isempty(par.type)
				error('Options ''order'' and ''type'' must be specified for method %s.',par.method)
			end
		case 'sysRnd'
			if isempty(par.scale) || isempty(par.order) || isempty(par.type)
				error('Options ''scale'', ''order'' and ''type'' must be specified for method %s.',par.method)
			end
		case 'scaleRnd'
			if isempty(par.scale)
				error('Option ''scale'' must be specified for method %s.',par.method)
			end
		otherwise
			error('Unsupported multipole type. Allowed are ''relToNom'', ''sysRnd'' or ''scaleRnd''.')
	end
end

% Applies the multipole to the lattice element
function RING = applyMultipoles(RING,ord,AB)
	if isfield(RING{ord},'PolynomAOffset')
		RING{ord}.PolynomAOffset = addPadded(RING{ord}.PolynomAOffset, AB(:,1)');
	else
		RING{ord}.PolynomAOffset = AB(:,1)';
	end

	if isfield(RING{ord},'PolynomBOffset')
		RING{ord}.PolynomBOffset = addPadded(RING{ord}.PolynomBOffset, AB(:,2)');
	else
		RING{ord}.PolynomBOffset = AB(:,2)';
	end

% 	RING{ord}.MaxOrder=length(RING{ord}.PolynomBOffset)-1;
end

% Auziliary function to add unequally long arrays (zero padding)
function v = addPadded(v1,v2)
	if ~((iscolumn(v1)&&iscolumn(v2))||(isrow(v1)&&isrow(v2)))
		error('Wrong dimensions.');
	end
	l1=length(v1);
	l2=length(v2);
	if l2>l1; v1(l2)=0; end
	if l1>l2; v2(l1)=0; end
	v=v1+v2;
end
