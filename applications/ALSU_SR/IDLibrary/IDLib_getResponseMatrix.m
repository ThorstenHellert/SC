function out = IDLib_getResponseMatrix(SC,quadOrds,evalOrds,varargin)
% IDLib_getResponseMatrix
% =======================
%
% NAME
% ----
% IDLib_getResponseMatrix - Calculates the respose matrix of the ID compensation figure of merrit
%
% SYNOPSIS
% --------
% `out = IDLib_getResponseMatrix(SC,quadOrds,evalOrds [,options])`
%
% DESCRIPTION
% -----------
% This function calcualtes the respose matrix of the ID compensation figure of merrit.
%
% INPUT
% -----
% `SC`::
% 	The SC base structure including not closed IDs and without errors.
% `quadOrds`::
% 	Quadrupole ordinates for which the response matrix should be calcualted.
% `evalOrds`::
% 	Ordinates at which the response matrix should be evaluated.
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'dK'` (`1E-2`)::
%	Quadrupole setpoint change.
% `'refOnly'` (`0`)::
%	If true, only the figure of merrit reference value is calcualted instead of the response matrix.
% `'muWeighting'` (`1`)::
%	Weighting factor for phase advance vs. beta function error.
% `'magnetWeighting'` (`1`)::
%	Weighting factors for individual magnets.
% `'vertWeighting'` (`1`)::
%	Weighting factors for vertical vs. horizontal plane.
%
% RETURN VALUE
% ------------
% `deltaK`::
% 	Vector of quadrupole setpoint changes.
% `finalFOM`::
% 	Value of final figure of merrit.
%
% SEE ALSO
% --------
% *IDlib_applyCorrection*, *IDlib_calcCorrection*


	% Parse optional arguments
	p = inputParser;
	addOptional(p,'dK',1E-2);
	addOptional(p,'refOnly',0);
	addOptional(p,'muWeighting',1);
	addOptional(p,'vertWeighting',1);
	addOptional(p,'magnetWeighting',1);
	parse(p,varargin{:});
	par = p.Results;

	
	% Get reference values
	[out,nu0,~] = atlinopt(SC.IDEALRING,0,evalOrds);
	beta0       = vertcat(out.beta);
	mu0         = vertcat(out.mu)/2/pi;

	% Get initial values
	[out,nu1,~] = atlinopt(SC.RING,0,evalOrds);
	beta1       = vertcat(out.beta);
	mu1         = vertcat(out.mu)/2/pi;

	% Get figure of merrit reference value
	ref = [((beta1-beta0)./beta0)                    .* repmat(par.magnetWeighting(:),1,2) .* repmat([1 par.vertWeighting],length(evalOrds),1) ;...
		   (mu1-mu0)              .* par.muWeighting .* repmat(par.magnetWeighting(:),1,2) .* repmat([1 par.vertWeighting],length(evalOrds),1);...
		    nu1-nu0 ];
	ref = ref(:);
	
	if par.refOnly
		out = ref;
	else 
		% Allocate matrix
		Mat = nan(4*length(evalOrds)+2,length(quadOrds));
		
		% Loop over quadrupoles
		for nQuad=1:length(quadOrds)
			
			% Change quadrupole value
			SC.RING{quadOrds(nQuad)}.PolynomB(2) = SC.RING{quadOrds(nQuad)}.PolynomB(2) + par.dK;
			
			% Get lattice values
			[out,nu,~] = atlinopt(SC.RING,0,evalOrds);			
			beta       = vertcat(out.beta);
			mu         = vertcat(out.mu)/2/pi;
	
			% Generate matrix
			tmp = ( [(beta1-beta)./beta1                     .* repmat(par.magnetWeighting(:),1,2) .* repmat([1 par.vertWeighting],length(evalOrds),1) ;...
					  (mu1-mu)            .* par.muWeighting .* repmat(par.magnetWeighting(:),1,2) .* repmat([1 par.vertWeighting],length(evalOrds),1);...
				       nu1-nu] ) / par.dK;
			
			Mat(:,nQuad) = tmp(:);
			
			% Re-set quadrupole value
			SC.RING{quadOrds(nQuad)}.PolynomB(2) = SC.RING{quadOrds(nQuad)}.PolynomB(2) - par.dK;
		end
	
		out = Mat;
	end
end
