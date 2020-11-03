function RING = SCgetModelRING(SC,varargin)
% SCgetModelRING
% ==============
%
% NAME
% ----
% SCgetModelRING - Returns a model lattice based on current setpoints
%
% SYNOPSIS
% --------
% `RING = SCgetModelRING(SC, [, options])`
%
%
% DESCRIPTION
% -----------
% This function calculates a model lattice based on the setpoints of `SC.RING`. Misalignments, 
% lattice errors and dipole fields are excluded.
%
%
% RETURN VALUES
% -------------
% `RING`::
%	The idealised RING structure
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'includeAperture'` (`0`)::
%	If true, the returned model ring includes the aperture 

	p = inputParser;
	addOptional(p,'includeAperture',0);
	parse(p,varargin{:});

	
	% Create local copy of ideal ring
	RING = SC.IDEALRING;
	
	for ord=1:length(SC.RING)
		% Get magnet setpoints
		if isfield(SC.RING{ord},'SetPointA') && isfield(SC.RING{ord},'SetPointB')
			RING{ord}.PolynomA = SC.RING{ord}.SetPointA;
			RING{ord}.PolynomB = SC.RING{ord}.SetPointB;
			
			% Remove dipole fields
			RING{ord}.PolynomA(1) = 0.0;
			RING{ord}.PolynomB(1) = 0.0;
		end
		
		% Check if aperture should be included
		if p.Results.includeAperture
			if isfield(SC.RING{ord},'EApertures')
				RING{ord}.EApertures = SC.RING{ord}.EApertures;
			end
			if isfield(SC.RING{ord},'RApertures')
				RING{ord}.RApertures = SC.RING{ord}.RApertures;
			end
		end

		% Get cavity setpoints
		if isfield(SC.ORD,'Cavity') && isfield(RING{ord},'Frequency')
			for field={'Frequency','Voltage','TimeLag'}
				RING{ord}.(field{1}) = SC.RING{ord}.([field{1} 'SetPoint']);
			end
		end
	end

end
