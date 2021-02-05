function RING = setRndPolynoms_ALSU_SR(RING)
% DESCRIPTION
% -----------
% Sets the random multipole errors for the ALS-U storage ring



	famNames = {'QD1','QD2','QD3','QF1','QF2','QF3','QF4','QF5','QF6'};
	for fam=famNames
		fprintf('Including RND multipoles for magnet: %s\n',fam{1})
		RING = applyScaledRandom(RING,fam{1},sprintf('%s_GAUSS_POLE_LENGTH_ASSEMBLY_PHY.tsv',fam{1}));
	end

	famNamesENG = {'SF','SD','SH1','SH2'};
	famNamesAT  = {'SF','SD','SHH$','SHH2'};
	for nFam=1:length(famNamesENG)
		fprintf('Including RND multipoles for magnet: %s\n',famNamesAT{nFam})
		RING = applyScaledRandom(RING,famNamesAT{nFam},sprintf('%s_GAUSS_POLE_LENGTH_ASSEMBLY_PHY.tsv',famNamesENG{nFam}));
	end

	fprintf('Including RND multipoles for magnet: HBEND\n')
	RING = applyScaledRandom(RING,'SB','SR_HBEND_RND_PHY.tsv');
	
	
	fprintf('Including RND multipoles for magnet: BEND1\n')
	RING = applyScaledRandom(RING,'BEND1','SR_BENDA_RND_AT-4Nov2020.tsv');
	
	famNames = {'BEND2','BEND3'};
	for fam=famNames
		fprintf('Including RND multipoles for magnet: %s\n',fam{1})
		RING = applyScaledRandom(RING,fam{1},'SR_BENDB_RND_AT-4Nov2020.tsv');
	end
	
end

function RING = applyScaledRandom(RING,re,fname)
	% Get ordinates of considered magnets
	ords=SCgetOrds(RING,re);
	% Read multipole errors from table
	AB = SCmultipolesRead(fname);
	
	fprintf('%d magnets found.\n',length(ords))
	
	% Account for bench measurements and correction/shimming
	if strncmp(re,'Q',1)
		fprintf('Dipole & quadrupole terms removed for %s.\n',fname)
		AB(1:2,:) = 0;
	elseif strncmp(re,'SB',2)
		fprintf('Dipole terms removed for %s.\n',fname)
		AB(1,:) = 0;
	elseif strncmp(re,'S',1)
		fprintf('Dipole & quadrupole & sextupole terms removed for %s.\n',fname)
		AB(1:3,:) = 0;
	end
	
	% Apply multipole errors to magnets
	RING = SCsetMultipoles(RING,ords,AB,'method','rnd');
end
