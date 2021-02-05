function RING = setRndPolynoms_ALSU_AR(RING,scale)
% Sets the random multipole errors for the ALS-U accumulator ring

	
	RING = applyScaledRandom(RING,':BEND:',scale,'AR_BEND_GAUSS_POLE_LENGTH_ASSEMBLY_GENERAL-SPEC_CASE#2_PHY.tsv');

	RING = applyScaledRandom(RING,':QD:',  scale,'AR_QD_GAUSS_POLE_LENGTH_ASSEMBLY_Jan_data_wd_4sig_cut_2sig_CLIP_PHY.tsv');
	RING = applyScaledRandom(RING,':QFA:', scale,'AR_QFA_GAUSS_POLE_LENGTH_ASSEMBLY_Jan_data_wd_4sig_cut_2sig_CLIP_PHY.tsv');
	RING = applyScaledRandom(RING,':QF:',  scale,'AR_QF_GAUSS_POLE_LENGTH_ASSEMBLY_Jan_data_wd_4sig_cut_2sig_CLIP_PHY.tsv');
	RING = applyScaledRandom(RING,':SD:',  scale,'AR_SD_GAUSS_POLE_LENGTH_ASSEMBLY_Jan_data_wd_4sig_cut_2sig_CLIP_PHY.tsv');
	RING = applyScaledRandom(RING,':SF:',  scale,'AR_SF_GAUSS_POLE_LENGTH_ASSEMBLY_Jan_data_wd_4sig_cut_2sig_CLIP_PHY.tsv');
	RING = applyScaledRandom(RING,':SHD:', scale,'AR_SHD_GAUSS_POLE_LENGTH_ASSEMBLY_Jan_data_wd_4sig_cut_2sig_CLIP_PHY.tsv');
	RING = applyScaledRandom(RING,':SHF:', scale,'AR_SHF_GAUSS_POLE_LENGTH_ASSEMBLY_Jan_data_wd_4sig_cut_2sig_CLIP_PHY.tsv');

end

function RING = applyScaledRandom(RING,re,scale,fname)
	% Get ordinates of considered magnets
	ords=SCgetOrds(RING,re);
	% Read multipole errors from table
	[AB,order,type] = SCmultipolesRead(fname);
		
	% Apply multipole errors to magnets
	RING = SCsetMultipoles_ALSU_AR(RING,ords,AB,'method','scaleRnd','order',order,'type',type,'scale',scale);
end
