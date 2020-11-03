function eta = SCgetModelDispersion(SC,BPMords,CAVords,varargin)
% SCgetModelDispersion
% ====================
%
% NAME
% ----
% SCgetModelDispersion - Calculates the lattice dispersion based on current setpints
%
% SYNOPSIS
% --------
% `eta = SCgetModelDispersion(SC,BPMords,CAVords [, options])`
%
%
% DESCRIPTION
% -----------
% Calcualtes the dispersion at the ordinates `BPMords` by changing the frequency of the rf cavities
% specified in `CAVords` using the current magnet setpoints without any roll/alignment/calibration
% errors. Optionally the design lattice is used.
%
%
% RETURN VALUES
% -------------
% `eta`::
%	The dispersion given in [m/Hz].
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'rfStep'` (1E3)::
%	Change of rf frequency [Hz]
% `'useIdealRing'` (0)::
%	If true, the design lattice specified in `SC.IDEALRING` is used.
%
% SEE ALSO
% --------
% *SCgetDispersion*, *SCgetModelRM*


	% Parse optional arguments
	p = inputParser;
	addOptional(p,'rfStep',1E3);
	addOptional(p,'useIdealRing',0);
	parse(p,varargin{:});
	par=p.Results;

	% Create local copy of ideal ring
	RING = SC.IDEALRING;

	% Check if only ideal ring should be used, otherwise set to Polynoms to the SetPoint (if it exists),
	% null the dipole components and remove offsets, rolls, apertures.
	if ~par.useIdealRing
		for ord=1:length(RING)
			if isfield(RING{ord},'SetPointA') && isfield(RING{ord},'SetPointB')
				% Read magnet setpoints
				RING{ord}.PolynomA = SC.RING{ord}.SetPointA;
				RING{ord}.PolynomB = SC.RING{ord}.SetPointB;

				% Remove dipole fields
				RING{ord}.PolynomA(1) = 0.0;
				RING{ord}.PolynomB(1) = 0.0;
			end
			% Remove offsets etc.
			rmfs = intersect(fieldnames(RING{ord}),{'T1','T2','R1','R2','EApertures','RApertures'});
			RING{ord} = rmfield(RING{ord},rmfs);
		end
	end

	% Number of BPMs and CMs
	nBPM = length(BPMords);

	% Allocate result matrix
	eta = nan(2 * nBPM , 1);

	% Calculate initial trajectories
	Bref = findorbit6(RING,BPMords);

	if any(isnan(Bref(:)))
		fprintf('Initial orbit is NaN. Aborting. \n')
		return
	end

	% Change RF frequency
	for ord = CAVords(:)'
		RING{ord}.Frequency = RING{ord}.Frequency + par.rfStep;
	end

	% Get second BPM reading
	B = findorbit6(RING,BPMords);

	% Calculate dispersion
	eta = reshape((B([1 3],:)'-Bref([1 3],:)')/par.rfStep,[],1);

end

