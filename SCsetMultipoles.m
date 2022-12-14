function RING = SCsetMultipoles(RING,ords,AB,varargin)
% SCsetMultipoles
% ===============
%
% NAME
% ----
% SCsetMultipoles - sets multipole errors in lattice elements
%
% SYNOPSIS
% --------
% `RING = SCsetMultipoles(RING, ords, AB [, options])`
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
% `AB`::    [N x 2] array of PolynomA/B multipole errors.
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'method'` (`'rnd'`)::
%   Specifies which multipole error method should be applied. Possible are
%   * `'sys'`:   This option is intended for normalized systematic multipole error tables. It sets 
%                the systematic multipoles of the field component defined by option `'order'`
%                and `'type'`. It is required that the `AB` entries are normalized by that component,
%                e.g. `AB(2,1)=1` for skew-quadrupole systematic multipoles. The systematic 
%                multipoles are from now on scaled with the current magnet excitation and added to the 
%                PolynomA/B fields.
%   * `'rnd'`:   This option is intended for random multipole error tables. It randomly generates
%                multipole components with a 2-sigma truncated Gaussian distribution from each of
%                the `AB` entries. The final multipole errors are stored in the PolynomA/BOffset of
%                the lattice elements.
% `'order'` ([])::
%   Numeric value defining the order of the considered magnet: [1,2,3,...] => [dip,quad,sext,...]
% `'type'` ([])::
%   Numeric value defining the type of the considered magnet: [1,2] => [skew/normal]
%
%
% EXAMPLES
% --------
% Defines random multipole components for the 'QF' magnet and adds it to the field offsets of all 
% magnets named 'QF'.
% ------------------------------------------------------------------
% ords = SCgetOrds(SC.RING,'QF');
% AB = [0 1E-5;...
%       0 1E-4;...
%       0 0;...
%       0 1E-2];
% RING = SCsetMultipoles(RING,ords,AB,'method','rnd');
% ------------------------------------------------------------------
%
% Reads the normalized systematic multipole components for the skew quadrupole excitation of the 'SF'
% magnet from table 'SF_skew_quad_AT_norm.tsv' and assigns it to all magnets named 'SF'. Note that 
% the data table in 'SF_skew_quad_AT_norm.tsv' must be 1.0 for the skew quadrupole component.
% ------------------------------------------------------------------
% ords = SCgetOrds(SC.RING,'SF');
% [AB,order,type] = SCmultipolesRead('SF_skew_quad_AT_norm.tsv');
% RING = SCsetMultipoles(RING,ords,AB,'method','sys','order',order,'type',type);
% ------------------------------------------------------------------
%
%
% SEE ALSO
% --------
% *SCmultipolesRead*, *SCupdateMagnets*



	p = inputParser;
	addParameter(p,'method','rnd');
	addOptional(p,'order',[]);
	addOptional(p,'type',[]);
	parse(p,varargin{:});
	par = p.Results;
	
	inputCheck(par);


	% Loop over magnets
	for ord=ords
		% Check which multipole error method should be applied
		switch par.method
			case 'sys'

				% Apply multipole errros to lattice
				RING = applySysMultipoles(RING,ord,AB,par.order,par.type);

			case 'rnd'

				% Generate random multipoles 
				addAB = SCrandnc(2,size(AB)) .* AB;

				% Apply multipole errros to lattice
				RING = applyRndMultipoles(RING,ord,addAB);
		end
	end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions


% Applies the systematic multipole to the lattice element
function RING = applySysMultipoles(RING,ord,AB,order,type)

	% Set nominal coefficent to zero (we dont want to add to it)
	AB(order,type) = 0;

	if type==1
		% Add systmatic PolynomA entries for coil A('order')
		RING{ord}.SysPolAFromA{order} = AB(:,1)';
		% Add systmatic PolynomB entries for coil A('order')
		RING{ord}.SysPolBFromA{order} = AB(:,2)';
	else
		% Add systmatic PolynomA entries for coil B('order')
		RING{ord}.SysPolAFromB{order} = AB(:,1)';
		% Add systmatic PolynomB entries for coil B('order')
		RING{ord}.SysPolBFromB{order} = AB(:,2)';
	end
	
end



% Applies the random multipole to the lattice element
function RING = applyRndMultipoles(RING,ord,AB)
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

end

% Auxiliary function to add unequally long arrays (zero padding)
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

% Check if input looks reasonable
function inputCheck(par)
	switch par.method
		case 'sys'
			if isempty(par.order) || isempty(par.type)
				error('Options ''order'' and ''type'' must be specified for method %s.',par.method)
			end
		case 'rnd'
			
		otherwise
			error('Unsupported multipole method. Allowed are ''sys'' or ''rnd''.')
	end
end
