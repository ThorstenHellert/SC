% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Workspace

clear all;clc;
%restoredefaultpath
global plotFunctionFlag SCinjections; 

% Initialize random number generator (important when running on cluster!)
devrand=fopen('/dev/urandom','r'); rngseed= fread(devrand,1,'*uint32'); rng(rngseed); fclose(devrand);
% Initialize results structure
results=struct('States',[]);runtime=tic;

% Set SC paths
mainDir = '~/sc';
addpath(mainDir);
% Set ALS-U AR paths
addpath(fullfile(mainDir,'applications/ALSU_AR_final/Multipoles'));
addpath(fullfile(mainDir,'applications/ALSU_AR_final/Lattices'));
addpath(fullfile(mainDir,'applications/ALSU_AR_final'));
% Set AT paths
addpath('~/at/atmat');
addpath('~/MML/applications/loco');
atpath()
% For convinience, switch off backtracing of warning messages
warning off backtrace


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Setup

nParticles = 100; % Number of particles in regular injections
nShots     = 1;   % Number of injections for averaging BPM readings

nParticlesQScan  = 100;  % Number of particles when performing tune scan
nRFParticles     = 100;  % Number of particles when performing calibrating RF cavity
nPartBeamCapture = 100;  % Number of particles when evaluating beam transmission
DAsteps          = 20;   % Number of angular steps for DA measurement

nTurnsBeamCapture = 500; % Number of turns when evaluating beam transmission
nTurnsTuneScan    = 500; % Number of turns when performing tune scan

postBBAoffset = 50E-6; % BPM offsets after pseudo-BBA
alpha         = 50;    % Regularization parameter for inverting response matrix
eps           = 1E-4;  % Noise threshold in trajectory/orbit feeback


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization

% Initialize toolkit
SC = SCinit(ALS_U_AR_v4_FINAL);

% Register AR
SC = register_ALSU_AR(SC);

% Define apertures
SC.RING = defineApertures_ALSU_AR(SC.RING);

% Get model response matrices
[RM1,RM2,MidealCO,Bref1,Bref2] = SCloadSaveMideal(SC,'AR');

% Apply errors
SC = SCapplyErrors(SC);

% Load systematic higher order multipole terms
SC.RING = setNominalPolynoms_ALSU_AR(SC.RING);
SC.RING = setSysPolynoms_ALSU_AR(SC.RING,100E-6,0.05);
SC.RING = setRndPolynoms_ALSU_AR(SC.RING,1);
SC = SCupdateMagnets(SC,SC.ORD.Magnet);

% Switch off cavity
SC.RING = SCcronoff(SC.RING,'cavityoff');

% Switch off sextupoles
sextOrds = SCgetOrds(SC.RING,'SF|SD|SHF|SHD');
SC = SCsetMags2SetPoints(SC,sextOrds,2,3,0,'method','abs');

% Set injection parameters
SC.INJ.nParticles = nParticles;
SC.INJ.nShots     = nShots;
SC.INJ.trackMode  = 'TBT';
SC.INJ.nTurns     = 1;

% Select if progress should be plotted
plotFunctionFlag = []; % Plot every injection
globPlotFlag     = 1; % Plot results at various stages


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Correction Chain

for FOO=0 % FAKE Loop. Just something to ``break'' out of if there is a fatal error.
		
	% Get pseudo inverse
	Minv1 = SCgetPinv(RM1,'alpha',alpha);
    
	results.States{end+1} = 'initial';
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% FIRST TURN 
	
	[CUR,ERROR] = SCfeedbackFirstTurn(SC,Minv1,'verbose',1);
	if ~ERROR; SC=CUR; results.States{end+1} = 'afterFirstTurn'; else; break; end
	
	[CUR,ERROR] = SCfeedbackRun(SC,Minv1,'target',300E-6,'maxsteps',30,'eps',eps,'verbose',1);
	if ~ERROR; SC=CUR; results.States{end+1} = 'afterFeedback'; else; break; end
	
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% GET 2-TURN TRANSMISSION
	SC.INJ.nTurns = 2;
	
	% Get pseudo inverse
	Minv2 = SCgetPinv(RM2,'alpha',alpha);
	
	for nBPM=[12 6 4 3 2]
		[CUR,ERROR] = SCfeedbackStitch( SC,Minv2,'nBPMs',nBPM,'maxsteps',100,'verbose',1);
		if ~ERROR; break; end
	end
	if ~ERROR
		SC=CUR; results.States{end+1} = 'afterStich'; 
	else
		B = SCgetBPMreading(CUR);
		lastBPM = find(isnan(B(1,length(SC.ORD.BPM)+1:end)),1);
		[CUR,ERROR] = SCfeedbackStitch(CUR,Minv2,'nBPMs',max(1,min(lastBPM-10,lastBPM)),'maxsteps',100,'verbose',1);

		if ~ERROR; SC=CUR; results.States{end+1} = 'afterStich'; else; break; end
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% MULTI-TURN ORBIT
	[CUR,ERROR] = SCfeedbackRun(SC,Minv2,'target',300E-6,'maxsteps',100,'eps',eps,'verbose',1);
	if ~ERROR; SC=CUR; results.States{end+1} = 'afterFeedback'; else; break; end
	
	[CUR,ERROR] = SCfeedbackBalance(SC,Minv2,'maxsteps',32,'eps',eps,'verbose',1);
	if ~ERROR; SC=CUR; results.States{end+1} = 'afterBalance'; else; break; end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Ramp up sextupoles 		
	nSteps=10;
	for nStep = linspace(1/nSteps,1,nSteps)
		% Increase sextupole strength
		SC = SCsetMags2SetPoints(SC,sextOrds,2,3,nStep,'method','rel');

		% Run feedback
		[CUR,ERROR] = SCfeedbackBalance(SC,Minv2,'maxsteps',10,'eps',eps,'verbose',1);
		if ~ERROR; SC=CUR; end
	end
	
	if ~ERROR; results.States{end+1} = 'afterSextupole'; else; break; end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Q-Scan to ensure transmission 
	
	qSPvec = {1 + linspace(-5E-2,5E-2,51) , 1 + linspace(-50E-2,50E-2,51)};
	qOrds  = SCgetOrds(SC.RING,{'QF:','QD'});
	[qSP,CUR,results.maxTurns,results.finTrans,ERROR] = SCtuneScan(SC,qOrds,qSPvec,'nParticles',nParticlesQScan,'nTurns',50,'target',0.8,'verbose',1,'plotFlag',globPlotFlag,'fullScan',0);
	if ERROR~=2; SC=CUR; results.States{end+1} = 'afterTransmission'; end

    
    SC.INJ.nParticles  = nRFParticles;
    SC.INJ.trackMode   = 'TBT';
    
    % Switch cavity on
    SC.RING = SCcronoff(SC.RING,'cavityon');
    
    % Store initial RF state
	tmpCO = findorbit6(SC.RING);
    results.rfEnergyError(1) = SC.INJ.Z0(5) - tmpCO(5);
    results.rfPhaseError(1)  = SC.INJ.Z0(6) - tmpCO(6);
	
	for nIter=1:3
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% RF phase correction 
		[deltaPhi,ERROR] = SCsynchPhaseCorrection(SC,...
												'nTurns',25,...
												'nSteps',25,...
												'plotResults',globPlotFlag,...
												'verbose',1);
		
		% Apply phase correction
		if ~ERROR; SC = SCsetCavs2SetPoints(SC,SC.ORD.Cavity,'TimeLag',deltaPhi,'add'); results.States{end+1} = 'afterPhase'; else; break; end
						
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% RF frequency correction 
		fMax = 1E3;
		[deltaF,ERROR] = SCsynchEnergyCorrection(SC,...
												'range',fMax*[-1 1],...
												'nTurns',130,...
												'minTurns',20,...
												'nSteps',15,...
												'plotResults',globPlotFlag,...
												'verbose',1);
		
		if ~ERROR
			% Check for boundaries (if phase correction went poorly, deltaF might be unreasonably large)
			if deltaF~=min(max(deltaF,-fMax),fMax)
				fprintf('Energy correction %.2fkHz above threshold and will be cliped.\n',fMax)
				deltaF = min(max(deltaF,-fMax),fMax);
			end
			% Apply frequency correction
			SC = SCsetCavs2SetPoints(SC,SC.ORD.Cavity,'Frequency',deltaF,'add');

			results.States{end+1} = 'afterFrequency';
		end
        
		 % Store current RF state
        tmpCO = findorbit6(SC.RING);
        results.rfEnergyError(1+nIter) = SC.INJ.Z0(5) - tmpCO(5);
        results.rfPhaseError(1+nIter)  = SC.INJ.Z0(6) - tmpCO(6);

    end
    
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Q-Scan to achieve beam capture %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	qSPvec = {1 + linspace(-5E-2,5E-2,51) , 1 + linspace(-50E-2,50E-2,51)};
	qOrds  = SCgetOrds(SC.RING,{'QF:','QD'});
	[qSP,CUR,results.maxTurns,results.finTrans,ERROR] = SCtuneScan(SC,qOrds,qSPvec,'nParticles',nParticlesQScan,'nTurns',nTurnsTuneScan,'target',0.8,'verbose',1,'plotFlag',globPlotFlag,'fullScan',0);
	if ERROR~=2; SC=CUR; results.States{end+1} = 'afterTransmission'; end

	results.qSPvec   = qSPvec;
	results.qSP      = qSP;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Beam capture achieved (?) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[results.finMaxTurns,...
	 results.finLostCount,~] = SCgetBeamTransmission(SC,'nParticles',nPartBeamCapture,'nTurns',nTurnsBeamCapture,'verbose',1,'plotFlag',globPlotFlag);
	
	if ERROR 
		fprintf('Beam capture not achieved (%d turns)\n',results.finMaxTurns); 
		if results.finMaxTurns < 50; fprintf('Insufficient transmission to continue.\n'); break; end
	end
	
	SC.INJ.trackMode = 'ORB';
	fprintf('Switch to orbit mode\n')
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Beam based alignment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	QuadOrds = getBPM2QuadPairing_ALSU_AR(SC);
	SC       = SCpseudoBBA(SC,SC.ORD.BPM,QuadOrds,postBBAoffset);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Orbit correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

	% Dispersion scaling
	scaleDisp = 1E8;
	% Get response matrix
	MCO = SCgetRespMat(SC,5E-5,SC.ORD.BPM,SC.ORD.CM);
	if any(isnan(MCO(:))) || any(isinf(MCO(:)))
		MCO = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'trackMode','ORB');
		if any(isnan(MCO(:)))
			MCO = SCgetModelRM(SC.IDEALRING,SC.ORD.BPM,SC.ORD.CM,'trackMode','ORB');
		end
	end
	% Get dispersion
	eta = SCgetDispersion(SC,1E3);
	if any(isnan(eta(:)))
		eta = SCgetModelDispersion(SC,SC.ORD.BPM,SC.ORD.Cavity);
		if any(isnan(eta(:)))
			eta = SCgetModelDispersion(SC,SC.ORD.BPM,SC.ORD.Cavity,'useIdealRing',1);
		end
	end	

	% Define regularization loop parameters
	alphaVec = [30 15 10:-1:1];

	fprintf('Initial closed orbit deviation (hor//ver):   %.0fum // %.0fum\n',1E6*sqrt(mean(SCgetBPMreading(SC).^2,2)))
	for	alpha=alphaVec
		fprintf('Orbit correction with alpha = %d\n',alpha)
				
		% Check if closed orbit can be found
		if ERROR==0 && ~any(isnan(findorbit6(SC.RING,1)))
			SC.INJ.trackMode = 'ORB';
			fprintf('Switch to orbit mode\n')
		else
			SC.INJ.trackMode  = 'pORB';
			SC.INJ.nTurns     = 100;
			SC.INJ.nParticles = nPartBeamCapture;
			fprintf('Switch to pseudo orbit mode\n')
		end

        % Get pseudo inverse
        MplusCO = SCgetPinv([MCO scaleDisp*eta],'alpha',alpha);
        
        % Run feedback
        [CUR,ERROR] = SCfeedbackRun(SC,MplusCO,'target',0,'maxsteps',50,'scaleDisp',scaleDisp,'verbose',1);

		if ~ERROR	
			% Check for sufficient beam transmission for orbit mode
			[maxTurns,~,ERROR] = SCgetBeamTransmission(CUR,'nParticles',nPartBeamCapture,'nTurns',50,'verbose',1,'plotFlag',globPlotFlag);
			
			% Check if transmission is below pseudo-orbit validity
			if maxTurns < 50
				fprintf('Transmission below 50 turns, abort orbit correction.\n')
				alpha = alphaVec(max([find(alphaVec==alpha)-1,1])); break
			end

			B0rms = sqrt(mean(SCgetBPMreading(SC).^2,2)); 
			Brms  = sqrt(mean(SCgetBPMreading(CUR ).^2,2)); 
			
			if mean(B0rms) < mean(Brms)
				fprintf('No further improvement with alpha = %d\n\n\n\n',alpha);
				alpha = alphaVec(max([find(alphaVec==alpha)-1,1])); break
			else
				SC = CUR;
				fprintf('CO mprovement with alpha = %d:\n hor: %.0fum -> %.0fum\n ver: %.0fum -> %.0fum\n',alpha,1E6*B0rms(1),1E6*Brms(1),1E6*B0rms(2),1E6*Brms(2))
			end
		else
			alpha = alphaVec(max([find(alphaVec==alpha)-1,1])); break
		end
	end
	fprintf('Final closed orbit deviation (hor//ver): %.1fum // %.1fum   with alpha = %d.\n',1E6*sqrt(mean(SCgetBPMreading(SC).^2,2)),alpha)
	warning('on','SC:CM1');
	results.States{end+1} = 'afterOrbitcorrection';
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get pre LOCO lattice properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
    results.preLOCO    = SCcalcLatticeProperties(SC,'DAsteps',DAsteps);
    results.preLOCO.SC = SC;
 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% LOCO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	try
		SC = performLOCO_ALSU_AR(SC);	
		results.States{end+1} = 'final';
    catch ME
		fprintf('LOCO did not work...\n')
		results.postLOCO = [];
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get post LOCO lattice properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	results.postLOCO    = SCcalcLatticeProperties(SC,'DAsteps',DAsteps);
    results.postLOCO.SC = SC;

     
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Perform RF-BEND correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	sigmaF = 1*1.1E3;
	
	targetF  = SC.IDEALRING{SC.ORD.Cavity}.Frequency + sigmaF * SCrandnc(2);
	results.targetF  = targetF;
	try
		[SC,out] = performBENDRFCorrection_ALSU_AR(SC,targetF,DAsteps);
		results.preBEND  = out.preBEND;
		results.postBEND = out.postBEND;
		results.States{end+1} = 'post-BEND-RF';
	catch ME
		fprintf('BEND did not work...\n')
		fprintf('%s/n',ME.message)
		results.preBEND  = [];
		results.postBEND = [];
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Record runtime and injections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	results.runtime     = toc(runtime);
	results.nInjections = SCinjections;

end

% If correction was aborted before 'preLOCO' state was recorded, save SC state.
if ~isfield(results,'preLOCO')
    results.breakSC = SC;
end

