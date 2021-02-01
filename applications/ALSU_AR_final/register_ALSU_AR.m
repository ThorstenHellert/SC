function SC = register_ALSU_AR(SC)
% Registers the ALS-U AR in the SC toolkit and defines uncertainties
	
	% Factors for switching off entire error types % % % % % % % % % % % % % % % % % %
	injErrorFactor  = 1;
	injJitterFactor = 1;
	magErrorFactor  = 1;
	diagErrorFactor = 1;
	RFerrorFactor   = 1;
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Add magnet errors to lattice 
	magOffset     = 50E-6  * [1 1 0] * magErrorFactor; % Offsets of magnets in x, y and z [m]
	magRoll       = 200E-6 * [1 0 0] * magErrorFactor; % Roll of magnets around z-, x- and y-axis [rad]
	magCal        = 1E-3             * magErrorFactor; % Relative magnet strength error
	
	SectionOffset = 100E-6 * [1 1 0] * magErrorFactor; % Offsets of sections in x, y and z [m]
	GirderOffset  = 50E-6  * [1 1 0] * magErrorFactor; % Offsets of girders in x, y and z [m]
	girderRoll    = 100E-6 * [1 0 0] * magErrorFactor; % Roll of girders around z-, x- and y-axis [rad]
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Add errors to diagnostic devices% 
	BPMcal        = 5E-2   * [1 1] * diagErrorFactor; % BPM calibration error 
	BPMoffset     = 500E-6 * [1 1] * diagErrorFactor; % BPM offset [m]
	BPMnoise      = 10E-6  * [1 1] * diagErrorFactor; % BPM noise [m]
	BPMnoiseCO    = 1E-6   * [1 1] * diagErrorFactor; % BPM noise for stored beam [m]
	BPMroll       = 4000E-6        * diagErrorFactor; % BPM roll around z-axis [rad]
	CMcal         = 5E-2           * diagErrorFactor; % CM calibration error
	
	CMlimit       = 0.2E-3 * [1 1]; % [HCM VCM] limit [rad]
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Circumfernece error 	
	circError     = 1e-06; % Circumeference error [rel]
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Add errors to RF cavity 
	RFfrequency     = 1E2  * RFerrorFactor;         % RF frequency error [Hz]
	RFvoltage       = 5E-3 * RFerrorFactor * 1.2E6; % RF voltage error [V]
	RftimeLag       = 1/4  * RFerrorFactor * 0.6;   % RF phase [m]
		
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Circ shift lattice to account for injection setup 
	SC.RING      = circshift(SC.RING     , 1-SCgetOrds(SC.RING     ,'^DK'));
	SC.IDEALRING = circshift(SC.IDEALRING, 1-SCgetOrds(SC.IDEALRING,'^DK'));
	% Switch on cavity and radiation
	SC.RING = SCcronoff(SC.RING,'cavityon','radiationon');
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Register SHF&SHD: HCM/VCM/Sext	
	ords = SCgetOrds(SC.RING,'SHF|SHD');
	SC = SCregisterMagnets(SC,ords,...
		'HCM',CMlimit(1),...
		'VCM',CMlimit(2),...
		'CalErrorB',[CMcal 0 magCal],...
		'CalErrorA',[CMcal 0 0     ],...
		'MagnetOffset',magOffset,...
		'MagnetRoll',magRoll);

	% Register SKEW SF: HCM/VCM/Skew/Sext	
	ords = SCgetOrds(SC.RING,'SF');
	ords = ords(3:4:end); % Odd sectors the first SF
	SC = SCregisterMagnets(SC,ords,...
		'HCM',CMlimit(1),...
		'VCM',CMlimit(2),...
		'SkewQuad',Inf,...
		'CalErrorB',[CMcal 0      magCal],...
		'CalErrorA',[CMcal magCal 0     ],...
		'MagnetOffset',magOffset,...
		'MagnetRoll',magRoll);

	% Register OTHER SF: HCM/VCM/Sext	
	ords = SCgetOrds(SC.RING,'SF');
	ords = setdiff(ords,ords(3:4:end));
	SC = SCregisterMagnets(SC,ords,...
		'HCM',CMlimit(1),...
		'VCM',CMlimit(2),...
		'CalErrorB',[CMcal 0 magCal],...
		'CalErrorA',[CMcal 0 0     ],...
		'MagnetOffset',magOffset,...
		'MagnetRoll',magRoll);

	% Register SKEW SD: Skew/Sext
	ords = SCgetOrds(SC.RING,'SD');
	ords = ords(1:4:end); % Even sectors the first SD
	SC = SCregisterMagnets(SC,ords,...
		'SkewQuad',Inf,...
		'CalErrorB',[0 0      magCal],...
		'CalErrorA',[0 magCal 0     ],...
		'MagnetOffset',magOffset,...
		'MagnetRoll',magRoll);

	% Register OTHER SD: Sext
	ords = SCgetOrds(SC.RING,'SD');
	ords = setdiff(ords,ords(1:4:end));
	SC = SCregisterMagnets(SC,ords,...
		'CalErrorB',[0 0 magCal],...
		'CalErrorA',[0 0 0     ],...
		'MagnetOffset',magOffset,...
		'MagnetRoll',magRoll);

	% Register BEND: Combined Function
	ords = SCgetOrds(SC.RING,'BEND');
	SC = SCregisterMagnets(SC,ords,...
		'CF',1); % No CalError because of multipole table / No offsets and rolls because ind. rafts

	% Register QF&QFA&QD: Quadrupoles
	ords = SCgetOrds(SC.RING,'QF:|QD|QFA');
	SC = SCregisterMagnets(SC,ords,...
		'CalErrorB',[0 magCal],...
		'MagnetOffset',magOffset,...
		'MagnetRoll',magRoll);

	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Register BPMs 
	ords = SCgetOrds(SC.RING,'BPM');
	SC = SCregisterBPMs(SC,ords,...
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
		'Offset',GirderOffset,...
		'Roll',girderRoll);

	% Register sections
	sectOrds = [SCgetOrds(SC.RING,'SECTIONSTART');SCgetOrds(SC.RING,'SECTIONEND')];
	SC = SCregisterSupport(SC,...
		'Section',sectOrds,...
		'Offset',SectionOffset);




	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Define initial bunch emittance
	emit0 = 300E-9*[1 0.1]; % [m*rad] -> 10% coupling
	% Define longitudinal size of initial bunch
	beamSize5 = 1E-3;   % Relative energy error
	beamSize6 = 2.5E-3; % Relative phase error [m]
	% Calculate Courant-Snyder parameters at beginning of ring
	TD = twissring(SC.IDEALRING,0,1);
	alpha = [TD.alpha(1) TD.alpha(2)];
	beta  = [TD.beta(1)  TD.beta(2)];
	% Define sigma matrices of injected beam
	sigx = emit0(1) * [beta(1) -alpha(1) ; -alpha(1) (1+alpha(1).^2)./beta(1)];
	sigy = emit0(2) * [beta(2) -alpha(2) ; -alpha(2) (1+alpha(2).^2)./beta(2)];
	sigz = [beamSize5 0 ; 0 beamSize6].^2;
	% Define 6x6 beam size matrix
	SC.INJ.beamSize = blkdiag(sigx,sigy,sigz);

	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Define random injection errors
	SC.SIG.randomInjectionZ(1) = injJitterFactor * 50E-6;
	SC.SIG.randomInjectionZ(2) = injJitterFactor * 5E-6; % 0.1% of fast kicker strength
	SC.SIG.randomInjectionZ(3) = injJitterFactor * 5E-6;
	SC.SIG.randomInjectionZ(4) = injJitterFactor * 2E-6;
	SC.SIG.randomInjectionZ(5) = injJitterFactor * 1E-4;
	SC.SIG.randomInjectionZ(6) = injJitterFactor * 0.1/360 * 299792458/SC.IDEALRING{SC.ORD.Cavity(1)}.Frequency;
	
	% Define systematic injection errors
	SC.SIG.staticInjectionZ(1) = injErrorFactor * 600E-6;
	SC.SIG.staticInjectionZ(2) = injErrorFactor * 150E-6;
	SC.SIG.staticInjectionZ(3) = injErrorFactor * 500E-6;
	SC.SIG.staticInjectionZ(4) = injErrorFactor * 100E-6;
	SC.SIG.staticInjectionZ(5) = injErrorFactor * 1E-3;
	SC.SIG.staticInjectionZ(6) = injErrorFactor * 0;
	
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Circumference uncertainty
	SC.SIG.Circumference = circError;
	
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Percentage of partcles which can be lost while still getting BPM reading
	SC.INJ.beamLostAt    = 0.6;

end
