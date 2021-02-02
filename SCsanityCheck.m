function SCsanityCheck(SC)
% SCsanityCheck
% =============
%
% NAME
% ----
% SCsanityCheck - Checks if the current registration looks reasonable
%
% SYNOPSIS
% --------
% `SCsanityCheck(SC)`
%
% DESCRIPTION
% -----------
% Performs a sanity check on the current `SC` structure adn returns warnings if things look fishy.
%
% INPUTS
% ------
% `SC`:: SC base structure
%
% SEE ALSO
% --------
% *SCregisterMagnets*, *SCregisterBPMs*, *SCregisterCAVs*


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Check if anything is registered
	if ~isfield(SC,'ORD')
		error('Nothing is registered.')
	else


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Check if BPMs are registered
	if ~isfield(SC.ORD,'BPM')
		warning('No BPMs registered. Use ''SCregisterBPMs''.')
	else
		% Check if BPMs are really registered
		if isempty(SC.ORD.BPM)
			warning('No BPMs registered. Use ''SCregisterBPMs''.');
		else
			fprintf('%d BPMs registered.\n',length(SC.ORD.BPM))
		end
		% Check if BPMs are uniquely registered
		if length(unique(SC.ORD.BPM))~=length(SC.ORD.BPM)
			warning('BPMs not uniquely defined.')
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Check if support structures are registered
	if ~isfield(SC.ORD,'Girder') && (isfield(SC.ORD,'Plinth') || isfield(SC.ORD,'Section'))
		warning('Girders must be registered for other support structure misalingments to work.')
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Check if CMs are registered
	if ~isfield(SC.ORD,'CM')
		warning('No CMs registered. Use ''SCregisterCMs''.')
	else
		if isempty(SC.ORD.CM{1})
			warning('No horizontal CMs registered. Use ''SCregisterCMs''.')
		else
			fprintf('%d HCMs registered.\n',length(SC.ORD.CM{1}))
		end
		if length(SC.ORD.CM)~=2 || isempty(SC.ORD.CM{2})
			warning('No vertical CMs registered. Use ''SCregisterCMs''.')
		else
			fprintf('%d VCMs registered.\n',length(SC.ORD.CM{2}))
		end
		% Check if CMs are uniquely registered
		if length(unique(SC.ORD.CM{1}))~=length(SC.ORD.CM{1})
			warning('Horizontal CMs not uniquely defined.')
		end
		if length(unique(SC.ORD.CM{2}))~=length(SC.ORD.CM{2})
			warning('Vertical CMs not uniquely defined.')
		end
		for ord=SC.ORD.CM{1}
			if SC.RING{ord}.CMlimit(1)==0
				warning('HCM limit is zero (Magnet ord: %d). Sure about that?',ord)
			end
		end
		for ord=SC.ORD.CM{2}
			if SC.RING{ord}.CMlimit(2)==0
				warning('VCM limit is zero (Magnet ord: %d). Sure about that?',ord)
			end
		end
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Check if magnets are registered correctly
	if ~isfield(SC.ORD,'Magnet')
		warning('No magnets are registered. Use ''SCregisterMagnets''.')
	else
		for ord=SC.ORD.Magnet
			% Check if lengths if normal/skew fields don't match
			if length(SC.RING{ord}.PolynomB)~=length(SC.RING{ord}.PolynomA)
				error('Length of PolynomB and PolynomA are not equal (Magnet ord: %d)',ord)
			elseif length(SC.RING{ord}.SetPointB)~=length(SC.RING{ord}.CalErrorB)
				warning('Length of SetPointB and CalErrorB are not equal (Magnet ord: %d)',ord)
			elseif length(SC.RING{ord}.SetPointA)~=length(SC.RING{ord}.CalErrorA)
				warning('Length of SetPointA and CalErrorA are not equal (Magnet ord: %d)',ord)
			end
			if isfield(SC.RING{ord},'PolynomBOffset')
				if length(SC.RING{ord}.PolynomBOffset)~=length(SC.RING{ord}.PolynomAOffset)
					error('Length of PolynomBOffset and PolynomAOffset are not equal (Magnet ord: %d)',ord)
				end
			end
			% Check if combined function looks reasonable
			if isfield(SC.RING{ord},'CombinedFunction') && SC.RING{ord}.CombinedFunction==1
				if ~isfield(SC.RING{ord},'BendingAngle')
					error('Combined function magnet (ord: %d) requires field ''BendingAngle''.',ord)
				end
				if SC.RING{ord}.NomPolynomB(2)==0 || SC.RING{ord}.BendingAngle==0
					warning('Combined function magnet (ord: %d) has zero bending angle or design quadrupole component.',ord)
				end
			end
			% Check if sigma structure has wrong fields
			if isfield(SC.SIG,'Mag') && ~isempty(SC.SIG.Mag{ord})
				for field=fieldnames(SC.SIG.Mag{ord})'
					if ~isfield(SC.RING{ord},field)
						warning('Field in SC.SIG.Mag doesn''t match lattice element (Magnet ord: %d)',ord)
					end
					if strcmp(field{1},'MagnetOffset')
						if length(SC.SIG.Mag{ord}.(field{1}))~=3
							warning('SC.SIG.Mag{%d}.MagnetOffset should be a [1x3] (dx,dy,dz) array.',ord)
						end
					end
				end
			end
			% Check if Master/Child definition looks reasonable
			if isfield(SC.RING{ord},'MasterOf')
				masterFields=fieldnames(SC.RING{ord});
				for cOrd=SC.RING{ord}.MasterOf
					for field=fieldnames(SC.RING{cOrd})'
						if ~any(strcmp(masterFields,field{1}))
							error('Child magnet (ord: %d) has different field ''%s'' than master magnet (ord: %d).',cOrd,field{1},ord)
						end
					end
				end
			end

		end
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Check if cavities are registered
	if ~isfield(SC.ORD,'Cavity')
		warning('No cavity registered. Use ''SCregisterCAVs''.')
	else
		% Check if cavities are really registered
		if isempty(SC.ORD.Cavity)
			warning('No cavity registered. Use ''SCregisterBPMs''.');
		else
			fprintf('%d cavity/cavities registered.\n',length(SC.ORD.Cavity))
		end
		% Check if cavities are uniquely registered
		if length(unique(SC.ORD.Cavity))~=length(SC.ORD.Cavity)
			warning('Cavities not uniquely defined.')
		end
		% Check if sigma structure has wrong fields
		if isfield(SC.SIG,'RF')
			for ord=SC.ORD.Cavity
				for field=fieldnames(SC.SIG.RF{ord})'
					if ~isfield(SC.RING{ord},field)
						warning('Field in SC.SIG.RF doesn''t match lattice element (Cavity ord: %d)',ord)
					end
				end
			end
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Check if beam properties are set correctly
	if any(size(SC.INJ.beamSize)~=[6,6])
		error('6x6 sigma matrix has to be used!')
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Check if apertures are set properly
	apEl=[];
	for ord=1:length(SC.RING)
		if isfield(SC.RING{ord},'EApertures') && isfield(SC.RING{ord},'RApertures')
			warning('Lattice element #%d has both EAperture and RAperture',ord)
		end
		if isfield(SC.RING{ord},'EApertures') || isfield(SC.RING{ord},'RApertures')
			apEl = [apEl ord];
		end
	end
	if isempty(apEl)
		fprintf('No apertures found.\n')
	else
		fprintf('Apertures defined in %d out of %d elements.\n',length(apEl),length(SC.RING))
	end

end
