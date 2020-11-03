function RING = setNominalPolynoms_ALSU_AR(RING)
% Sets the systematic multipole fields for the ALS-U accumulator ring


	% Primary coils
	RING = applyRelative(RING,':BEND:','AR_BEND_systematic_PHY_NORM.tsv');
	RING = applyRelative(RING,':QD:',  'AR_QD_systematic_PHY_NORM.tsv');
	RING = applyRelative(RING,':QF:',  'AR_QF_systematic_PHY_NORM.tsv');
	RING = applyRelative(RING,':QFA:', 'AR_QFA_systematic_PHY_NORM.tsv');
	RING = applyRelative(RING,':SD:',  'AR_SD_systematic_sext_PHY_NORM.tsv');
	RING = applyRelative(RING,':SF:',  'AR_SF_systematic_sext_PHY_NORM.tsv');
	RING = applyRelative(RING,':SHD:', 'AR_SHD_systematic_sext_PHY_NORM.tsv');
	RING = applyRelative(RING,':SHF:', 'AR_SHF_systematic_sext_PHY_NORM.tsv');
		
end

function RING = applyRelative(RING,re,fname)
	% Get ordinates of considered magnets
	ords=SCgetOrds(RING,re);
	% Read multipole errors from table
	[AB,order,type] = SCmultipolesRead(fname);
	% Set dipole component to zero
	if order~=1; AB(1,:)=0; end
	% Apply multipole errors to magnets
	RING = SCsetMultipoles(RING,ords,AB,'method','relToNom','order',order,'type',type);
end
