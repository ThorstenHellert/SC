function RING = setSysPolynoms_ALSU_SR(RING)
	% Sets the systematic multipoles for the ALSU-SR	

	
	
	% Nominal multipoles
	RING = applySys(RING,'QD1'    ,'QD1_SYS_fromAL1344-8997_PHY_NORM.tsv');
	RING = applySys(RING,'QD2'    ,'QD3_SYS_fromAL1344-8997_PHY_NORM.tsv'); % QD3 magnet tables for QD2 magnets
	RING = applySys(RING,'QD3'    ,'QD3_SYS_fromAL1344-8997_PHY_NORM.tsv');
	RING = applySys(RING,'QF1'    ,'QF1_SYS_fromAL1344-8997_PHY_NORM.tsv');
	RING = applySys(RING,'QF2'    ,'QF2_SYS_fromAL1344-8997_PHY_NORM.tsv');
	RING = applySys(RING,'QF3'    ,'QF3_SYS_fromAL1344-8997_PHY_NORM.tsv');
	RING = applySys(RING,'QF4'    ,'QF4_SYS_fromAL1344-8997_PHY_NORM.tsv'); 
	RING = applySys(RING,'QF5'    ,'QF4_SYS_fromAL1344-8997_PHY_NORM.tsv'); % QF4 magnet tables for QF5 magnets
	RING = applySys(RING,'QF6'    ,'QF4_SYS_fromAL1344-8997_PHY_NORM.tsv'); % QF4 magnet tables for QF6 magnets
	RING = applySys(RING,'BEND1'  ,'BEND1_SYS_fromAL1344-8997_PHY_NORM.tsv');
	RING = applySys(RING,'BEND2'  ,'BEND2_SYS_fromAL1344-8997_PHY_NORM.tsv');
	RING = applySys(RING,'BEND3'  ,'BEND3_SYS_fromAL1344-8997_PHY_NORM.tsv');
	RING = applySys(RING,'SB'     ,'HBEND_SYS_fromAL1344-8997_PHY_NORM.tsv');
	RING = applySys(RING,'SF'     ,'SF_SYS_fromAL1344-8997_PHY_NORM.tsv');
	RING = applySys(RING,'SD'     ,'SD_SYS_fromAL1344-8997_PHY_NORM.tsv');
	RING = applySys(RING,'SHH$'   ,'SHH_SYS_fromAL1344-8997_PHY_NORM.tsv');
	RING = applySys(RING,'SHH2'   ,'SHH2_SYS_fromAL1344-8997_PHY_NORM.tsv');

	% Dipole correctors coils
	RING = applySys(RING,'QF2'    ,'QF2_Bx_corrector_130AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'QF2'    ,'QF2_By_corrector_130AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'QF3'    ,'QF3_Bx_corrector_215AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'QF3'    ,'QF3_By_corrector_215AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'QF[4-6]','QFRC_Bx_corrector_85AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'QF[4-6]','QFRC_By_corrector_85AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	
	RING = applySys(RING,'BEND2','B2_dipole_corrector_1795AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'BEND3','B3_dipole_corrector_1927AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');

	RING = applySys(RING,'SD','SF_Bx_corrector_50AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv'); % SF magnet tables for SD magnets
	RING = applySys(RING,'SD','SF_By_corrector_55AT_27p5AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'SF'  ,'SF_Bx_corrector_50AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'SF'  ,'SF_By_corrector_55AT_27p5AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'SHH$','SH1_Bx_corrector_1190AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'SHH$','SH1_By_corrector_3090AT_1545AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'SHH2','SH2_Bx_corrector_1190AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'SHH2','SH2_By_corrector_1378AT_689AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	
	% Skew-quad correctors coils
	RING = applySys(RING,'SF'  ,'SF_skew_quad_corrector_240AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'SHH2','SH2_skew_quad_corrector_336AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
	RING = applySys(RING,'SHH2','SH2_skew_quad_corrector_336AT_effective_multipoles_by_physics_unit_PHY_NORM.tsv');
end

function RING = applySys(RING,re,fname)
	% Get ordinates of considered magnets
	ords=SCgetOrds(RING,re);

	% Read multipole errors from table
	[AB,order,type] = SCmultipolesRead(fname);
	
	% BEND2/3 HCM coil will be calibrated with main coil such that no additional quadrupole component is excited
	if ~isempty(regexp(fname,'B2_dipole_corrector')) || ~isempty(regexp(fname,'B3_dipole_corrector'))
		fprintf('Zeroing quadrupole terms for HCM in %s!\n',re)
		AB(2,:) = 0;
	end
	
	% Apply multipole errors to magnets
	RING = SCsetMultipoles(RING,ords,AB,'method','sys','order',order,'type',type);
end
