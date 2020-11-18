function SC = SCpseudoBBA(SC,BPMords,MagOrds,postBBAoffset,varargin)
% SCpseudoBBA
% ===========
%
% NAME
% ----
% SCpseudoBBA - mimics the results of a BBA procedure
%
% SYNOPSIS
% --------
% `[SC] = SCpseudoBBA(SC, BPMords, MagOrds, postBBAoffset)`
%
%
% DESCRIPTION
% -----------
% Sets the BPM offsets to values as would be achieved by performing a BBA procedure. The BPMs 
% defined in `BPMords` are paired with the magnets defined in `MagOrds` and the BPM offsets 
% are adjusted to match the magnet offsets with an accuracy of `postBBAoffset`.
%
%
% INPUT
% -----
% `SC`::
%	SC base structure.
% `BPMords`::
%	[2 x N] array of BPM ordinates.
% `MagOrds`::
%	[2 x N] array of magnet ordinates.
% `postBBAoffset`::
%	[2 x N] array or single value of rms accuracy of pseudo BBA procedure [m].
%
%
% OPTIONS
% -------
% The following options can be specified as name-value pairs:
%
% `'sigma'` (`2`)::
%	Cutoff (number of sigmas) for generating the final offset distribution.
%
%
% RETURN VALUES
% -------------
% `SC`::
%	SC-structure with adjusted BPM offsets.
%
%
% EXAMPLES
% --------
%
% Assuming that each QF, QD and QFA magnet is accompanied by one BPM, the following assumes a 
% sucessful BBA routine and assigns a 50um BPM offset. 
% ------------------------------------------------------------------
% QuadOrds = repmat(SCgetOrds(SC.RING,'QF|QD|QFA'),2,1);
% SC = SCpseudoBBA(SC,SC.ORD.BPM,QuadOrds,50E-6);
% ------------------------------------------------------------------
%
% SEE ALSO
% --------
%  *SCgetOrds*, *SCBBA*



% Parse optional arguments
	p = inputParser;
	addOptional(p,'sigma',2);
	parse(p,varargin{:});

	if isempty(postBBAoffset)
		fprintf('No BPM offset defined. NO BBA PERFORMED!.\n')
		return
	else
		fprintf('Performing pseudo BBA with stored beam. New BPM offset is %.0fum (rms) w.r.t. quadrupole centers.\n',1E6*sqrt(mean(postBBAoffset(:).^2)))
	end
	
	if length(postBBAoffset)==1
		postBBAoffset = repmat(postBBAoffset,2,size(BPMords,2));
	end
	
	
	
	% Loop over BPMs
	for nBPM=1:size(BPMords,2)

		% Horizontal/vertical
		for nDim=1:2
		
			% Write new BPM offset		
			SC.RING{BPMords(nDim,nBPM)}.Offset(nDim) = SC.RING{MagOrds(nDim,nBPM)}.MagnetOffset(nDim) + SC.RING{MagOrds(nDim,nBPM)}.SupportOffset(nDim) - SC.RING{BPMords(nDim,nBPM)}.SupportOffset(nDim) + postBBAoffset(nDim,nBPM) * SCrandnc(p.Results.sigma);
		
		end
	end
	
	
end
