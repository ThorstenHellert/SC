% Clear memory and comand line window
clear all;clc;
% Restore default matlab path to delete any paths from other scripts
restoredefaultpath;
% Set the plotFunctionFlag and number of injections as global variables
global plotFunctionFlag SCinjections runParallel;
results=struct('States',[]);runtime=tic;

% Start with a random seed (on PC not really an issue, on the cluster however it happens that all the seeds start with the same random number generator)
devrand=fopen('/dev/urandom','r'); rngseed= fread(devrand,1,'*uint32'); rng(rngseed); fclose(devrand);

% Define the main SC directory (the one where you put all the downloaded files)
addpath('~/sc/');

% Define main SR directory 
mainDir = '~/sc/applications/ALSU_SR';
% Include the folder where SR routines and scripts (e.g. SR error definitions) are stored
addpath(mainDir);
% Include the folder where all SR lattice files are stored
addpath(fullfile(mainDir,'lattices'));
addpath(fullfile(mainDir,'Multipoles/SR_SYS_20210112_converted/'));
addpath(fullfile(mainDir,'Multipoles/SR_RND_20200825_converted/'));
addpath(fullfile(mainDir,'Multipoles/SR_CM_SYS_200903_converted/'));
addpath(fullfile(mainDir,'Multipoles/SR_RND_HBEND_converted/'));
addpath(fullfile(mainDir,'Multipoles/SR_RND_BEND_converted/'));
% Include ID library
addpath(fullfile(mainDir,'IDLibrary'));
addpath(fullfile(mainDir,'IDLibrary/kickmaps_COSMIC'));
addpath(fullfile(mainDir,'IDLibrary/kickmaps_LEDA'));
addpath(fullfile(mainDir,'IDLibrary/kickmaps_TENDER'));
addpath(fullfile(mainDir,'IDLibrary/kickmaps_XType'));
addpath(fullfile(mainDir,'IDLibrary/kickmaps_EPU35'));
addpath(fullfile(mainDir,'IDLibrary/kickmaps_EPU90'));
addpath(fullfile(mainDir,'IDLibrary/kickmaps_EPU50'));

% Put the folder where the function atpath.m is located to the matlab path
addpath('~/at/atmat');
% Put the folder where the function loco.m is located to the matlab path
addpath('~/MML/applications/loco/');
% Let AT set its paths
atpath()

% For convinience, switch off back-tracing warning messages
warning off backtrace

% Set IDs
results.IDs        = { 'COSMIC'          , 'LEDA'          , 'TENDER'           , 'U114' , 'XType'            ,'EPU35'               ,'EPU90'                 ,'EPU50'                };
results.IDsect     = [   5               ,   9             ,    1               ,  10    ,    6               ,  4                   ,  3                     ,  11                   ];
results.IDmap      = { 'CPShFCShFC.xls'  , 'Kickg4p2.xls'  ,  'TenderKick.xls'  ,  ''    ,  'XtypeIPKick.xls' ,'EPU35-shim6b-HR.prn' , 'EPU90-shim4b-g22.prn' , 'EPU50-good-shims-v1.prn' };

% Load lattice
RING = ALS_U_v21_4raft_SB;

% Remove first BPM after injection (space constraints)
tmp = SCgetOrds(RING,'BPM1$');
RING(tmp(1)) = [];
% Remove last BPM before injection (space constraints)
tmp = SCgetOrds(RING,'BPM19');
RING(tmp(end)) = [];

% Include IDs
RING = IDlib_includeIDs(RING,results.IDs,'Section',results.IDsect);

% Initialize toolkit
SC = SCinit(RING);

% Register ALSU-SR
[SC,BPMords,CMords] = register_ALSU_SR(SC);
% Save ideal SC state for ID compensation calculation
results.SCrefID = SC;
% Save BPM and CM ords used in orbit correction
results.BPMords = BPMords;
results.CMords  = CMords;

% Define apertures
SC.RING = setApertures_ALSU_SR(SC.RING);

% Get 1- and 2-turn model response matrices
RM1 = SCgetModelRM(SC,BPMords,CMords,'nTurns',1,'useIdealRing',1);
RM2 = SCgetModelRM(SC,BPMords,CMords,'nTurns',2,'useIdealRing',1);

% Apply errors
SC = SCapplyErrors(SC);

% Load multipole errors
SC.RING = setSysPolynoms_ALSU_SR(SC.RING);
SC.RING = setRndPolynoms_ALSU_SR(SC.RING);
SC = SCupdateMagnets(SC);

% Set static magnet field strength error
for ord=setdiff(SC.ORD.Magnet,SCgetOrds(SC.RING,'SB'))
	ind = find(SC.RING{ord}.NomPolynomB);
	if isfield(SC.RING{ord},'PolynomBOffset')
		SC.RING{ord}.PolynomBOffset(ind) = SC.RING{ord}.PolynomBOffset(ind) + 1E-3 * SCrandnc(2) * SC.RING{ord}.NomPolynomB(ind);
	else
		SC.RING{ord}.PolynomBOffset(ind) = 1E-3 * SCrandnc(2) * SC.RING{ord}.NomPolynomB(ind);
		SC.RING{ord}.PolynomAOffset      = zeros(size(SC.RING{ord}.PolynomBOffset));
	end
end
SC = SCupdateMagnets(SC);

% Define quadrupole trim coil vs sextupole center error
for ord = SCgetOrds(SC.RING,'SF|SD')
	SC.RING{ord}.QuadSextOffsetError = 6E-6 * SCrandnc(2,1,2);
end

% Switch off cavity
SC.RING = SCcronoff(SC.RING,'cavityoff');

% Switch off sextupoles
sextOrds = SCgetOrds(SC.RING,'SF|SD|SHH$|SHH2');
SC = SCsetMags2SetPoints(SC,sextOrds,2,3,0,'method','abs');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following parameters allow for convinient switching between 'debug' mode and 'real' calculation

runParallel        = 0; % Run parfor (in some functions)
plotFunctionFlag   = 0; % Switch on/off plotting every beam
globPlotFlag       = 0; % Switch on/off plotting of function results

% Switch on/off time consuming things
performRealCOBBA   = 1;  % Switch on/off stored beam BBA
performRealtrajBBA = 1;  % Switch on/off trajectory BBA
debugMode          = 1;  % Switch on/off debug mode
DAsteps            = 50; % Number of steps for DA calculation

% LOCO parameters
LOCOSVmethod = 1E-2;
LOCOcmStep   = 1 * 50E-6;
LOCOrfStep   = 1 * 3e2 * LOCOcmStep/1E-5;
RMpar        = {'mode','fixedKick'  ,'nSteps',3,'verbose',1}; % For response matrix measurement
LOCOhorDispWeight = 5;
LOCOverDispWeight = 5;

% Fake BBA values
postBBAoffsetTraj(1) = 1* 150E-6;
postBBAoffsetTraj(2) = 1* 50E-6;
postBBAoffsetQuad    = 1* 20E-6;
postBBAoffsetSext    = 1* 15E-6;


if debugMode
	nPartBeamCapture  = 5;  % Turns for final transmission
	nTurnsBeamCapture = 10; % Number of particles to check for beam capture
	SC.INJ.nParticles = 1;  % Particles per shot
	nPartQScan        = 10; % Number of particles to perform tune scan
else
	nTurnsBeamCapture = 1000; % Turns for final transmission
	nPartBeamCapture  = 50;   % Number of particles to check for beam capture
	SC.INJ.nParticles = 200;  % Particles per shot
	nPartQScan        = 200;  % Number of particles to perform tune scan
end

for FOO=0 % FAKE Loop. Just something to ``break'' out of if there is a fatal error.
	
	% Define regularization parameter
	alpha1 = 100;
	alpha2 = 10;
	
	% Define noise level
	eps = 1E-4;
		
	% Get pseudo inverse
	Minv1 = SCgetPinv(RM1,'alpha',alpha1);
	
	% Write succesfull correction steps in structure to see where it failed
	results.States{end+1} = 'initial';
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% FIRST TURN 
	[CUR,ERROR] = SCfeedbackFirstTurn(SC,Minv1,'wiggleAfter',10,'verbose',1,'BPMords',BPMords,'CMords',CMords,'wiggleSteps',16);
	if ~ERROR; SC=CUR; results.States{end+1} = 'afterFirstTurn'; else; break; end
	
	[CUR,ERROR] = SCfeedbackRun(SC,Minv1,'target',300E-6,'maxsteps',30,'eps',eps,'verbose',1,'BPMords',BPMords,'CMords',CMords);
	if ~ERROR; SC=CUR; results.States{end+1} = 'afterFeedback'; else; break; end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COOL MULTI TURN ORBIT 
	
	SC.INJ.nTurns = 2;
	
	% Get pseudo inverse
	Minv2 = SCgetPinv(RM2,'alpha',alpha1);
	
	for nBPM=[15 10 5 3]
		[CUR,ERROR] = SCfeedbackStitch( SC,Minv2,'nBPMs',nBPM,'maxsteps',200,'verbose',1,'BPMords',BPMords,'CMords',CMords);
		if ~ERROR; break; end
	end
	if ~ERROR; SC=CUR; results.States{end+1} = 'afterStich'; else; break; end
	
	[CUR,ERROR] = SCfeedbackRun(SC,Minv2,'target',300E-6,'maxsteps',100,'eps',eps,'verbose',1,'BPMords',BPMords,'CMords',CMords);
	if ~ERROR; SC=CUR; results.States{end+1} = 'afterFeedback'; else; break; end

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% >>> Evaluate beam transmission <<<
    [results.trans.preTraBBA.turns,...
	 results.trans.preTraBBA.count,~] = SCgetBeamTransmission(SC,...
														'nParticles',nPartBeamCapture,...
														'nTurns',nTurnsBeamCapture,...
														'verbose',1,...
														'plotFlag',globPlotFlag);

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BBA 
	QuadOrds = getBPM2QuadPairing_ALSU_SR(SC,repmat(SC.ORD.BPM,2,1));
	% Get quadrupole setpoint variation
	magSPvec = repmat({linspace(0.95,1.05,2)},size(QuadOrds));
	magSPvec(:,ismember(QuadOrds(1,:),SCgetOrds(SC.RING,'QD1|QF1')	)) = {linspace(0.65,1.05,2)};

	% BPM offsets
	results.BPMOffsets{1} = cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'Offset'))' + cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'SupportOffset'))' - cell2mat(atgetfieldvalues(SC.RING(getBPM2QuadPairing_ALSU_SR(SC,SC.ORD.BPM)),'MagnetOffset',{1,1:2}))' - cell2mat(atgetfieldvalues(SC.RING(getBPM2QuadPairing_ALSU_SR(SC,SC.ORD.BPM)),'SupportOffset',{1,1:2}))';

	
	% Trajectory BBA and injection error correction loop
	for nIter=1:2

		if performRealtrajBBA
			qOrd = SCgetOrds(SC.RING,'QF');
			[SC,errorFlags] = SCBBA(SC,repmat(SC.ORD.BPM,2,1),QuadOrds,...
									'mode','TBT',...
									'fakeMeasForFailures',1,...
									'outlierRejectionAt',Inf,...
									'quadOrdPhaseAdvance',qOrd(1),...
									'quadStrengthPhaseAdvance',[0.95 0.8 1.05],...
									'magSPvec',magSPvec,...
									'nSteps',10,...
									'verbose',1,...
									'plotResults',globPlotFlag,...
									'plotLines',0,...
									'minSlopeForFit',0.03,...
									'minBPMrangeAtBBABBPM',500E-6,...
									'BBABPMtarget',1E-3,...
									'minBPMrangeOtherBPM',100E-6,...
									'maxStdForFittedCenters',600E-6);

		else
			fprintf('Performing pseudo-BBA.\n')
			SC = SCpseudoBBA(SC,repmat(SC.ORD.BPM,2,1),QuadOrds,postBBAoffsetTraj(nIter));
		end

		% BPM offsets
		results.BPMOffsets{end+1} = cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'Offset'))' + cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'SupportOffset'))' - cell2mat(atgetfieldvalues(SC.RING(getBPM2QuadPairing_ALSU_SR(SC,SC.ORD.BPM)),'MagnetOffset',{1,1:2}))' - cell2mat(atgetfieldvalues(SC.RING(getBPM2QuadPairing_ALSU_SR(SC,SC.ORD.BPM)),'SupportOffset',{1,1:2}))';




		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% GET CLOSED ORBIT BUMP

		% Get pseudo inverse
		Minv2 = SCgetPinv(RM2,'alpha',alpha2);

		[CUR,ERROR] = SCfeedbackRun(SC,Minv2,'target',50E-6,'maxsteps',100,'eps',eps,'verbose',1,'BPMords',BPMords,'CMords',CMords);
		if ~ERROR; SC=CUR; results.States{end+1} = 'afterFeedback'; else; break; end

		[CUR,ERROR] = SCfeedbackBalance(SC,Minv2,'maxsteps',32,'eps',eps,'verbose',1,'BPMords',BPMords,'CMords',CMords);
		if ~ERROR; SC=CUR; results.States{end+1} = 'afterBalance'; else; break; end



		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% >>> Evaluate beam transmission <<<
		[results.trans.(sprintf('postTraBBA%d',nIter)).turns,...
		 results.trans.(sprintf('postTraBBA%d',nIter)).count,~] = SCgetBeamTransmission(SC,...
																						'nParticles',nPartBeamCapture,...
																						'nTurns',nTurnsBeamCapture,...
																						'verbose',1,...
																						'plotFlag',globPlotFlag);




		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% INJECTION CORRECTION 

		% Save initial state
		results.INJ.(sprintf('preInjCorr%d',nIter)) = SC.INJ.Z0;

		[deltaZ0,ERROR] = SCfitInjectionZ(SC,'injectionDrift','verbose',1);

		if ~ERROR
			initialZ0 = SC.INJ.Z0;
			for nStep=1:5
				SC.INJ.Z0 = initialZ0 + deltaZ0 * nStep/5;
				[CUR,ERROR] = SCfeedbackRun(SC,Minv2,'target',50E-6,'maxsteps',30,'eps',eps,'verbose',1,'BPMords',BPMords,'CMords',CMords);
				if ~ERROR; SC=CUR; else; break; end
			end
			results.States{end+1} = 'afterInjection';
		else; break; end

		% Save results
		results.INJ.(sprintf('postInjCorr%d',nIter)) = SC.INJ.Z0;


		% Run feedback
		[SC,~] = SCfeedbackRun(SC,Minv2,'target',0,'maxsteps',100,'eps',eps,'verbose',1,'BPMords',BPMords,'CMords',CMords);

		
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% >>> Evaluate beam transmission <<<
    [results.trans.postInj.turns,...
	 results.trans.postInj.count,~] = SCgetBeamTransmission(SC,...
														'nParticles',nPartBeamCapture,...
														'nTurns',nTurnsBeamCapture,...
														'verbose',1,...
														'plotFlag',globPlotFlag);

	
	
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Ramp up sextupoles 
		
	nSteps=10;
	for nStep = linspace(1/nSteps,1,nSteps)
		% Increase sextupole strength
		SC = SCsetMags2SetPoints(SC,sextOrds,2,3,nStep,'method','rel');

		% Run feedback
		[CUR,ERROR] = SCfeedbackBalance(SC,Minv2,'maxsteps',10,'eps',eps,'verbose',1,'BPMords',BPMords,'CMords',CMords);
		if ~ERROR; SC=CUR; end
		
		
		% >>> Evaluate beam transmission <<<
		[results.trans.sextRamp(round(nStep*nSteps)).turns,...
		results.trans.sextRamp(round(nStep*nSteps)).count,~] = SCgetBeamTransmission(SC,...
																'nParticles',nPartBeamCapture,...
																'nTurns',nTurnsBeamCapture,...
																'verbose',1,...
																'plotFlag',globPlotFlag);
	end
	
	results.States{end+1} = 'afterSextupole';
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% >>> Evaluate beam transmission <<<
    [results.trans.postSext.turns,...
	 results.trans.postSext.count,~] = SCgetBeamTransmission(SC,...
														'nParticles',nPartBeamCapture,...
														'nTurns',nTurnsBeamCapture,...
														'verbose',1,...
														'plotFlag',globPlotFlag);



	% Switch cavity on
	SC.RING = SCcronoff(SC.RING,'cavityon');
	
	% Save initial energy error
	tmpCO = findorbit6(SC.RING);
    results.rfEnergyError(1) = SC.INJ.Z0(5) - tmpCO(5);
    results.rfPhaseError(1)  = SC.INJ.Z0(6) - tmpCO(6);

	for nIter=1:3			
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% RF phase correction 
		[deltaPhi,ERROR] = SCsynchPhaseCorrection(SC,...
												'nTurns',35,...
												'nSteps',25,...
												'plotResults',globPlotFlag,...
												'plotProgress',0,...
												'verbose',1);
		
		% Apply phase correction
		if ~ERROR; SC = SCsetCavs2SetPoints(SC,SC.ORD.Cavity,'TimeLag',deltaPhi,'add'); results.States{end+1} = 'afterPhase'; end
		if globPlotFlag;SCplotPhaseSpace(SC,'nParticles',10,'nTurns',300);end
				
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% RF frequency correction 
		if nIter==3
			fMax = 1E3;
			nTRF = 200;
			nMin = 120;
		else
			fMax = 1E3;
			nTRF = 100;
			nMin = 80;
		end
		
		if nIter==1
			% Define relative quadrupole setpoints
			qSPvec = {1 + linspace(-1E-2,1E-2,51) , 1 + linspace(-1E-2,1E-2,51)};
			% Define quadrupole ordinates
			qOrds  = SCgetOrds(SC.RING,{'^QF1','^QD1'});
			% Perform tune scan
			[qSP,CUR,maxTurns,finTrans,ERROR] = SCtuneScan(SC,qOrds,qSPvec,...
															'nParticles',nPartQScan,...
															'nTurns',100,...
															'target',0.8,...
															'verbose',1,...
															'plotFlag',globPlotFlag);
			if ERROR~=2; SC=CUR; results.States{end+1} = 'afterTuneScan'; end
		end
		
		
		[deltaF,ERROR] = SCsynchEnergyCorrection(SC,...
												'range',fMax*[-1 1],...
												'nTurns',nTRF,...
												'minTurns',nMin,...
												'nSteps',15,...
												'plotResults',globPlotFlag,...
												'plotProgress',globPlotFlag,...
												'verbose',1);
		
		if ~ERROR
			% Check for boundaries (if phase correction went poorly, deltaF might be unreasonably large)
			if deltaF~=min(max(deltaF,-fMax),fMax)
				fprintf('Energy correction %.0fkHz above threshold and will be cliped.\n',fMax)
				deltaF = min(max(deltaF,-fMax),fMax);
			end
			% Apply frequency correction
			SC = SCsetCavs2SetPoints(SC,SC.ORD.Cavity,'Frequency',deltaF,'add');
			results.States{end+1} = 'afterFrequency';
		end
		
		if globPlotFlag;SCplotPhaseSpace(SC,'nParticles',10,'nTurns',300);end
		
		% Save current longitudinal injection error
		tmpCO = findorbit6(SC.RING);
        results.rfEnergyError(1+nIter) = SC.INJ.Z0(5) - tmpCO(5);
        results.rfPhaseError(1+nIter)  = SC.INJ.Z0(6) - tmpCO(6);
	end
	

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get post-RF lattice properties
	
	 % Get post-RF lattice properties and save SC
	 results.postRF    = calcLatticeProperties_ALSU_SR(SC,'DAsteps',DAsteps);
     results.postRF.SC = SC;
  
	 % >>> Evaluate beam transmission <<<
	 [results.trans.postRF.turns,...
	  results.trans.postRF.count,~] = SCgetBeamTransmission(SC,...
														'nParticles',nPartBeamCapture,...
														'nTurns',nTurnsBeamCapture,...
														'verbose',1,...
														'plotFlag',globPlotFlag);

	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Switch to orbit mode and calculate responsa matrix and dispersion
	
	fprintf('Switch to orbit mode\n')
	SC.INJ.trackMode = 'ORB';
	
	% Save all relevant information for orbit feedback (e.g. in BBA)
	RM                  = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'trackMode','ORB','useIdealRing',1);	
	% Take for OFB relevant entries of response matrix
	RMstruct.RM         = RM([find(ismember(SC.ORD.BPM,BPMords)) length(SC.ORD.BPM)+find(ismember(SC.ORD.BPM,BPMords))],[find(ismember(SC.ORD.CM{1},CMords{1})) length(SC.ORD.CM{1})+find(ismember(SC.ORD.CM{2},CMords{2}))]);
	RMstruct.eta        = SCgetModelDispersion(SC,BPMords,SC.ORD.Cavity,'rfStep',5);
	RMstruct.scaleDisp  = 1E7;
	RMstruct.BPMords    = BPMords;
	RMstruct.CMords     = CMords;
	% Pseudo-inverse for BBA orbit bump
	RMstruct.MinvCO = SCgetPinv([RMstruct.RM RMstruct.scaleDisp*RMstruct.eta],'alpha',5);
	% Regularization parameters for orbit feedback loops
	RMstruct.alphaVec = [100 50 20 10 5:-1:1 0.9:-0.1:0];

	% Get CM steps for LOCO RM measurement (mean BPM amplitude of 50um)
	LOCOcmSteps{1} = 50E-6 ./ mean(abs(RM(1:length(SC.ORD.BPM),1:length(SC.ORD.CM{1}))));
	LOCOcmSteps{2} = 50E-6 ./ mean(abs(RM(length(SC.ORD.BPM)+1:end,length(SC.ORD.CM{1})+1:end)));



	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Beam based alignment
	
	% Sextupole BBA BPMs 
	BPMordsSextBBA  = repmat(SCgetOrds(SC.RING,'BPM4|BPM5|BPM15|BPM16'),2,1);
	% Quadrupole BBA BPMs 
	BPMordsQuadBBA  = repmat(setdiff(BPMords,BPMordsSextBBA(1,:)),2,1);
	% BBA quadrupoles
	QuadOrds = getBPM2QuadPairing_ALSU_SR(SC,BPMordsQuadBBA);
	% BBA sextupoles
	sextOrds    = repmat(SCgetOrds(SC.RING,'SF|SD'),2,1);
	% Get quadrupole setpoint variation
	magSPvec = repmat({linspace(0.95,1.05,2)},size(BPMordsQuadBBA));
	magSPvec(:,ismember(QuadOrds(1,:),SCgetOrds(SC.RING,'QD1')	)) = {linspace(0.8,1.05,2)};
	% Get quadrupole trim coil setpoint variation
	magSPvecSext = repmat({linspace(-0.3,0.3,2)},size(BPMordsSextBBA));
	% Save BPM offset w.r.t. sextupole centers
	results.BPMOffsetsSext{1} = cell2mat(atgetfieldvalues(SC.RING(BPMordsSextBBA(1,:)),'Offset'))' + cell2mat(atgetfieldvalues(SC.RING(BPMordsSextBBA(1,:)),'SupportOffset'))'  - cell2mat(atgetfieldvalues(SC.RING(sextOrds(1,:)),'MagnetOffset',{1,1:2}))' - cell2mat(atgetfieldvalues(SC.RING(sextOrds(1,:)),'SupportOffset',{1,1:2}))';



	% Number of BBA/OFB/LOCO iterations
	for nIter=1:3
			
		if performRealCOBBA

			% Orbit correction (without weights)
			SC = performOrbitCorr_ALSU_SR(SC,RMstruct,'weight',[]);
			
			if nIter==1
				% >>> Evaluate beam transmission <<<
				[results.trans.preBBA.turns,...
				 results.trans.preBBA.count,~] = SCgetBeamTransmission(SC,...
																	'nParticles',nPartBeamCapture,...
																	'nTurns',nTurnsBeamCapture,...
																	'verbose',1,...
																	'plotFlag',globPlotFlag);										
			end
										
			% BBA on quadrupoles
			SC = SCBBA(SC,BPMordsQuadBBA,QuadOrds,...
							'mode','ORB',...
							'fakeMeasForFailures',1,...
							'outlierRejectionAt',200E-6,...
							'RMstruct',RMstruct,...
							'plotLines',0,...
							'plotResults',globPlotFlag,...
							'verbose',1,...
							'nSteps',10,...
							'magSPvec',magSPvec,...
							'minSlopeForFit',0.005,...
							'minBPMrangeAtBBABBPM',1*20E-6,...
							'BBABPMtarget',5*50E-6,...
							'minBPMrangeOtherBPM',0);
			

			% Orbit correction (without weights)
			SC = performOrbitCorr_ALSU_SR(SC,RMstruct,'weight',[]);

									
			% BBA on sextupoles (quad trim coils)					
			SC = SCBBA(SC,BPMordsSextBBA,sextOrds,...
							'mode','ORB',...
							'fakeMeasForFailures',1,...
							'outlierRejectionAt',200E-6,...
							'RMstruct',RMstruct,...
							'plotLines',0,...
							'plotResults',globPlotFlag,...
							'verbose',1,...
							'nSteps',10,...
							'switchOffSext',1,...
							'magSPflag','abs',...
							'magSPvec',magSPvecSext,...
							'minSlopeForFit',0.005,...
							'minBPMrangeAtBBABBPM',1*10E-6,...
							'BBABPMtarget',5*50E-6,...
							'minBPMrangeOtherBPM',0);
						
						
			% Add systematic quad trim-coil / sextupole center error
			for nBPM = 1:size(BPMordsSextBBA,2)
				SC.RING{BPMordsSextBBA(1,nBPM)}.Offset = SC.RING{BPMordsSextBBA(1,nBPM)}.Offset + SC.RING{sextOrds(1,nBPM)}.QuadSextOffsetError;
			end	
		else
			
			% BBA on quadrupoles
			SC = SCpseudoBBA(SC,BPMordsQuadBBA,QuadOrds,postBBAoffsetQuad);
			
			% BBA on sextupoles (quad trim coils)
			SC = SCpseudoBBA(SC,BPMordsSextBBA,sextOrds,postBBAoffsetSext);
		end
		
		% Orbit correction (with weights)
		SC = performOrbitCorr_ALSU_SR(SC,RMstruct);

		if nIter==1
			% >>> Evaluate beam transmission <<<
			[results.trans.preLOCO.turns,...
			 results.trans.preLOCO.count,~] = SCgetBeamTransmission(SC,...
																'nParticles',nPartBeamCapture,...
																'nTurns',nTurnsBeamCapture,...
																'verbose',1,...
																'plotFlag',globPlotFlag);
		end

		% Save lattices for orbit plot
		results.preLOCO{nIter}.RING = SC.RING;

		% Optics correction
		SC = performLOCO_ALSU_SR(SC,BPMords,CMords,...
									'RMpar',RMpar,...
									'diagCorr',1,...
									'CMstep',LOCOcmSteps,...
									'deltaRF',LOCOrfStep,...
									'RMstruct',RMstruct,...
									'HorizontalDispersionWeight',LOCOhorDispWeight,...
									'VerticalDispersionWeight'  ,LOCOverDispWeight);


		% Save BPM offsets w.r.t. quad centers
		results.BPMOffsets{end+1} = cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'Offset'))' + cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'SupportOffset'))' - cell2mat(atgetfieldvalues(SC.RING(getBPM2QuadPairing_ALSU_SR(SC,SC.ORD.BPM)),'MagnetOffset',{1,1:2}))' - cell2mat(atgetfieldvalues(SC.RING(getBPM2QuadPairing_ALSU_SR(SC,SC.ORD.BPM)),'SupportOffset',{1,1:2}))';
		% Save BPM offset w.r.t. sextupole centers
		results.BPMOffsetsSext{end+1} = cell2mat(atgetfieldvalues(SC.RING(BPMordsSextBBA(1,:)),'Offset'))' + cell2mat(atgetfieldvalues(SC.RING(BPMordsSextBBA(1,:)),'SupportOffset'))'  - cell2mat(atgetfieldvalues(SC.RING(sextOrds(1,:)),'MagnetOffset',{1,1:2}))' - cell2mat(atgetfieldvalues(SC.RING(sextOrds(1,:)),'SupportOffset',{1,1:2}))';

		results.tmpBBASC{nIter} = SC;
	end
	

	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get post LOCO lattice properties 

	% >>> Evaluate beam transmission <<<
    [results.trans.postLOCO.turns,...
	 results.trans.postLOCO.count,~] = SCgetBeamTransmission(SC,...
														'nParticles',nPartBeamCapture,...
														'nTurns',nTurnsBeamCapture,...
														'verbose',1,...
														'plotFlag',globPlotFlag);
	 
	 % Get post-LOCO lattice properties and save SC
	 results.postLOCO    = calcLatticeProperties_ALSU_SR(SC,'DAsteps',DAsteps);
     results.postLOCO.SC = SC;	
     
     % Save postLOCO RING for MA calcualtion
	 RING = SC.RING;
	 save('RING_postLOCO.mat','RING')
	 
	 
	 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% ID Compensation 
	try
		% Calculate ID correction
		[deltaK,finalFOM] = IDlib_calcCorrection(results.SCrefID,results.IDs,...
			'IDmap',results.IDmap,...
			'verbose',1);

		% Apply ID correction
		SC = IDlib_applyCorrection(SC,results.IDs,...
			'IDmap',results.IDmap,...
			'deltaK',deltaK,...
			'RMstruct',RMstruct,...
			'verbose',1);
				
	catch ME
		fprintf('ID compensartion did not work...\n')
		fprintf('%s/n',ME.message)
		results.postID = [];
		break;
	end
	
	% Get post-ID lattice properties and save SC
	results.postID    = calcLatticeProperties_ALSU_SR(SC,'DAsteps',DAsteps);
	results.postID.SC = SC;

	% Save postID RING for MA calcualtion
	RING = SC.RING;
	save('RING_postID.mat','RING')

end

% Save runtime and injections
results.runtime = toc(runtime);
results.nInjections = SCinjections;

if ~isfield(results,'postLOCO')
    results.breakSC = SC;
end
