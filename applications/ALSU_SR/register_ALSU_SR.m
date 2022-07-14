function [SC, BPMords, CMords] = register_ALSU_SR(SC)
	
	% Factors for switching off entire error types % % % % % % % % % % % % % % % % % %
	injErrorFactor  = 1;
	injJitterFactor = 1;
	magErrorFactor  = 1;
	diagErrorFactor = 1;
	RFerrorFactor   = 1;
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Add magnet errors to lattice 
	magOffset     = 1* [35E-6 35E-6 200E-6]   * magErrorFactor; % Offsets of magnets in x, y and z [m]
	magRoll       = 1* [200E-6 200E-6 200E-6] * magErrorFactor; % Roll of magnets around z-, x- and y-axis [rad]
	magCal        = 1* 2E-4                   * magErrorFactor; % Relative magnet strength error

	SectionOffset = 1* 100E-6 * [1 1 0]       * magErrorFactor; % Offsets of sections in x, y and z [m]
	GirderOffset  = 1* 35E-6  * [1 1 0;1 1 0] * magErrorFactor; % Offsets of girders in x, y and z [m]
	girderRoll    = 1* [100E-6 0 0]           * magErrorFactor; % Roll of girders around z-, x- and y-axis [rad]

	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Add errors to diagnostic devices% 
	BPMcal          = 1* 5E-2   * [1 1] * diagErrorFactor; % BPM calibration error 
	BPMoffset       = 1* 500E-6 * [1 1] * diagErrorFactor; % BPM offset [m]
	BPMnoise        = 1* 30E-6  * [1 1] * diagErrorFactor; % BPM noise [m]
	BPMnoiseCO      = 1* 1E-6   * [1 1] * diagErrorFactor; % BPM noise for stored beam [m]
	BPMroll         = 1* 4E-3           * diagErrorFactor; % BPM roll around z-axis [rad]
	BPMsumError     = 1* 0.2            * diagErrorFactor; % BPM sum signal calibration error
	CMcal           = 1* 5E-2           * diagErrorFactor; % CM calibration error
	
	CMlimitQuad     = 1* 0.8E-3;
	CMlimitSext     = 1* 0.4E-3;
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Circumfernece error 	
	circError       = 1* 1E-6; % Circumeference error [rel]
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Add errors to RF cavity 
	RFfrequency     = 1* 1E2  * RFerrorFactor; % RF frequency error [Hz]
	RFvoltage       = 1* 1E-3 * RFerrorFactor * SC.IDEALRING{findcells(SC.RING,'Frequency')}.Voltage; % RF voltage error [V]
	RftimeLag       = 1* 1/4  * RFerrorFactor * 299792458/SC.IDEALRING{findcells(SC.RING,'Frequency')}.Frequency; % RF phase [m]
	
	% Relative error of injected beam alpha and beta
	twissError      = 0.2; 
	

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get CMs and BPMs used for orbit correction
	ords = SCgetOrds(SC.RING,'^QF3');
	tmp  = SCgetOrds(SC.RING,'^QF4|^QF5|^QF6');ords = sort([ords tmp]);
	tmp  = SCgetOrds(SC.RING,'SHH$|SD|SHH2');ords = sort([ords tmp]);
	tmp  = SCgetOrds(SC.RING,'SF');ords = sort([ords tmp]);
	
	CMords  = {ords,ords};
	BPMords = setdiff(SCgetOrds(SC.RING,'BPM'),SCgetOrds(SC.RING,'BPM7$|BPM12$'));
	
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Register SB: Superbends 
	ords = SCgetOrds(SC.RING,'^SB$');
	SC = SCregisterMagnets(SC,ords,...
		'MagnetOffset',{magOffset,3},...
		'MagnetRoll',magRoll); 

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Register QD1/2/3b: Pure quad
	ords = SCgetOrds(SC.RING,'QD2b|QD3b');
	SC = SCregisterMagnets(SC,ords,...
		'CalErrorB',[0 magCal],...
		'MagnetOffset',{magOffset,3},...
		'MagnetRoll',magRoll); 
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Register QF1/QD1: Quad 
	ords = SCgetOrds(SC.RING,'^QF1|^QD1');
	SC = SCregisterMagnets(SC,ords,...
		'CalErrorB',[0 magCal],...
		'CalErrorA',[0 0     ],...
		'MagnetOffset',{magOffset,3},...
		'MagnetRoll',magRoll);


	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Register QF2/QF3: HCM/VCM/Quad (CF)
	ords = SCgetOrds(SC.RING,'^QF2|^QF3');
	SC = SCregisterMagnets(SC,ords,...
		'HCM',CMlimitQuad,...
		'VCM',CMlimitQuad,...
		'CF',1,...
		'CalErrorB',[CMcal magCal],...
		'CalErrorA',[CMcal 0     ],...
		'MagnetOffset',{magOffset,3},...
		'MagnetRoll',magRoll);

	% Register QF4-6: HCM/VCM/Quad (CF)
	ords = SCgetOrds(SC.RING,'^QF4|^QF5|^QF6');
	SC = SCregisterMagnets(SC,ords,...
		'HCM',CMlimitQuad,...
		'VCM',CMlimitQuad,...
		'CF',1,...
		'CalErrorB',[CMcal magCal],...
		'CalErrorA',[CMcal 0     ],...
		'MagnetOffset',{magOffset,3},...
		'MagnetRoll',magRoll);


	% Register SF/SD/SHF/SDF: HCM/VCM/Skew/Sext	
	ords = SCgetOrds(SC.RING,'SF|SD');
	SC = SCregisterMagnets(SC,ords,...
		'HCM',CMlimitSext,...
		'VCM',CMlimitSext,...
		'SkewQuad',Inf,...
		'CalErrorB',[CMcal 0      magCal],...
		'CalErrorA',[CMcal magCal 0     ],...
		'MagnetOffset',{magOffset,3},...
		'MagnetRoll',magRoll);

	% Register SF/SD/SHF/SDF: HCM/VCM/Skew/Sext	
	ords = SCgetOrds(SC.RING,'SHH$|SHH2');
	SC = SCregisterMagnets(SC,ords,...
		'HCM',CMlimitSext,...
		'VCM',CMlimitSext,...
		'SkewQuad',Inf,...
		'CalErrorB',[CMcal 0      magCal],...
		'CalErrorA',[CMcal magCal 0     ],...
		'MagnetOffset',{magOffset,3},...
		'MagnetRoll',magRoll);
	
	
	% Register BEND1: Combined Function
	ords = SCgetOrds(SC.RING,'BEND1');
	SC = SCregisterMagnets(SC,ords,...
		'CF',1,...
		'CalErrorB',[0 magCal],...
		'MagnetOffset',{magOffset,3},...
		'MagnetRoll',magRoll); 

	% Register BEND2-3: HCM/Combined Function
	ords = SCgetOrds(SC.RING,'BEND2|BEND3');
	SC = SCregisterMagnets(SC,ords,...
		'HCM',3.2E-3,...
		'CF',1,...
		'CalErrorB',[0 magCal],...
		'MagnetOffset',{magOffset,3},...
		'MagnetRoll',magRoll); 

	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Register BPMs 
	ords = SCgetOrds(SC.RING,'BPM');	
	SC = SCregisterBPMs(SC,ords,...
		'SumError',BPMsumError,...
		'CalError',BPMcal,...
		'Offset',BPMoffset,...
		'Roll',BPMroll,...
		'Noise',BPMnoise,...
		'NoiseCO',BPMnoiseCO);
		
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Register cavities
	CAVords = findcells(SC.RING,'Frequency');
	SC = SCregisterCAVs(SC,CAVords,...
		'FrequencyOffset',RFfrequency,...
		'VoltageOffset',RFvoltage,...
		'TimeLagOffset',RftimeLag);


	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Register girders
	girderOrds = reshape(SCgetOrds(SC.RING,'GIRDER'),2,[]);
	SC = SCregisterSupport(SC,...
		'Girder',girderOrds,...
		'Offset',{GirderOffset,3},...
		'Roll',girderRoll);


	% Register sections
	sectOrds              = SCgetOrds(SC.RING,'SECT');
 	sectOrds(2,:)         = circshift(sectOrds-1,-1);
	sectOrds(sectOrds==0) = length(SC.RING);	
	SC = SCregisterSupport(SC,...
		'Section',sectOrds,...
		'Offset',SectionOffset);

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Switch on cavity and radiation
	SC.RING = SCcronoff(SC.RING,'cavityon','radiationon');

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Define initial bunch emittance
	emit0 = 2E-9*[1 0.01]; % [m*rad] -> 1% coupling
	% Define longitudinal size of initial bunch
	beamSize5 = 1E-3;				% Relative energy error
	beamSize6 = 15E-12 * 299792458; % Relative phase error [m]
	% Calculate Courant-Snyder parameters at beginning of ring
	TD = twissring(SC.IDEALRING,0,1);
	alpha = [TD.alpha(1) TD.alpha(2)] .* (1 + twissError * SCrandnc(2,1,2));
	beta  = [TD.beta(1)  TD.beta(2)]  .* (1 + twissError * SCrandnc(2,1,2));

	% Define sigma matrices of injected beam
	sigx = emit0(1) * [beta(1) -alpha(1) ; -alpha(1) (1+alpha(1).^2)./beta(1)];
	sigy = emit0(2) * [beta(2) -alpha(2) ; -alpha(2) (1+alpha(2).^2)./beta(2)];
	sigz = [beamSize5 0 ; 0 beamSize6].^2;
	% Define 6x6 beam size matrix
	SC.INJ.beamSize = blkdiag(sigx,sigy,sigz);

	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Define random injection errors
	SC.SIG.randomInjectionZ(1) = injJitterFactor * 10E-6;
	SC.SIG.randomInjectionZ(2) = injJitterFactor * 10E-6; 
	SC.SIG.randomInjectionZ(3) = injJitterFactor * 1E-6;
	SC.SIG.randomInjectionZ(4) = injJitterFactor * 0.5E-6;
	SC.SIG.randomInjectionZ(5) = injJitterFactor * 1E-4;
	SC.SIG.randomInjectionZ(6) = injJitterFactor * 0.1/360 * 299792458/SC.IDEALRING{SC.ORD.Cavity(1)}.Frequency;
	
	% Define systematic injection errors
	SC.SIG.staticInjectionZ(1) = injErrorFactor * 500E-6;
	SC.SIG.staticInjectionZ(2) = injErrorFactor * 200E-6;
	SC.SIG.staticInjectionZ(3) = injErrorFactor * 500E-6;
	SC.SIG.staticInjectionZ(4) = injErrorFactor * 200E-6;
	SC.SIG.staticInjectionZ(5) = injErrorFactor * 1E-3;
	SC.SIG.staticInjectionZ(6) = injErrorFactor * 0;
	
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Circumference uncertainty
	SC.SIG.Circumference = circError;
	
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Percentage of partcles which can be lost while still getting BPM reading
	SC.INJ.beamLostAt    = 0.6;

end
