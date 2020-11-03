function RING = setSysPolynoms_ALSU_AR(RING,sigma_corr,sigma_skew)
% Sets the pseude random multipole fields induced by dipole and skew quadrupole correctors for the 
% ALS-U accumulator ring.

		
	% Calculate kick in PolynomA/B field value from length of CMs
	sigma_corr = sigma_corr/2.0300e-01;
	
	% Dipole correctors coils
	RING = applyRandom(RING,':SD:',  sigma_corr,'AR_SD_systematic_Bx_PHY_NORM.tsv');
	RING = applyRandom(RING,':SD:',  sigma_corr,'AR_SD_systematic_By_PHY_NORM.tsv');
	RING = applyRandom(RING,':SF:',  sigma_corr,'AR_SF_systematic_Bx_PHY_NORM.tsv');
	RING = applyRandom(RING,':SF:',  sigma_corr,'AR_SF_systematic_By_PHY_NORM.tsv');
	RING = applyRandom(RING,':SHD:', sigma_corr,'AR_SHD_systematic_Bx_PHY_NORM.tsv');
	RING = applyRandom(RING,':SHD:', sigma_corr,'AR_SHD_systematic_By_PHY_NORM.tsv');
	RING = applyRandom(RING,':SHF:', sigma_corr,'AR_SHF_systematic_Bx_PHY_NORM.tsv');
	RING = applyRandom(RING,':SHF:', sigma_corr,'AR_SHF_systematic_By_PHY_NORM.tsv');
	
	% Skew-quad correctors coils
	RING = applyRandom(RING,':SD:',  sigma_skew,'AR_SD_systematic_skewquad_PHY_NORM.tsv');
	RING = applyRandom(RING,':SF:',  sigma_skew,'AR_SF_systematic_skewquad_PHY_NORM.tsv');
	
end

function RING = applyRandom(RING,re,sigma,fname)
	% Get ordinates of considered magnets
	ords=SCgetOrds(RING,re);
	% Read multipole errors from table
	[AB,order,type] = SCmultipolesRead(fname);
	% Apply multipole errors to magnets
	RING = SCsetMultipoles(RING,ords,AB,'method','sysRnd','order',order,'type',type,'scale',sigma);
end
