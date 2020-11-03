function SC = performLOCO_ALSU_AR(SC)
% Performs the linear optics correction for the ALS-U accumulator ring
	
	
	%% Setup LOCO data structures
	
	CMstep   = 5E-5; % CM step for response matrix calculation
	deltaRF  = 1e2 * CMstep/1E-5; % RF difference for dispersion measurement
	
	% Setup LOCO model
	[LOCOmodel,LocoFlags,Init] = SClocoLib('setupLOCOmodel',SC,...
		'HorizontalDispersionWeight',1E1,...
		'VerticalDispersionWeight',2E1,...
		'SVmethod',1E-3);

	% Set up BPM and CM data structure structures
	[BPMData,CMData] = SClocoLib('getBPMCMstructure',SC,CMstep);
	
	% Get orbit response matrix and dispersion from measurement
	LocoMeasData = SClocoLib('getMeasurement',SC,CMstep,deltaRF,SC.ORD.BPM,SC.ORD.CM);
	
	% Initialize correction progress structure (for debugging)
	CorrProgress{1} = SC.RING;
	
	%% Run LOCO ++ first loop
	for nIter=1:7
		
		if nIter==1
			% Set up LOCO fit parameter structure for QF/QD magnets
			FitParameters = SClocoLib('setupFitparameters',SC,Init.SC.RING,LOCOmodel,deltaRF,...
				{SCgetOrds(SC.RING,'QF:'),'normal','individual',2E-3},... % {Ordinates, normal/skew, individual/family, deltaK}
				{SCgetOrds(SC.RING,'QD:'),'normal','individual',1E-3});   % {Ordinates, normal/skew, individual/family, deltaK}
		elseif nIter==2
			% Fit calibration factors from now on
			BPMData.FitGains    = 'Yes';
			CMData.FitKicks     = 'Yes';
		elseif nIter==4
			% Include off-diagonal ORM elements and dispersion, exclude calibration factors
			LocoFlags.Coupling   = 'Yes'; 
			LocoFlags.Dispersion = 'Yes'; 
			BPMData.FitGains     = 'No';
			CMData.FitKicks      = 'No';
			
			% Set up LOCO fit parameter structure for all magnets
			FitParameters = SClocoLib('setupFitparameters',SC,Init.SC.RING,LOCOmodel,deltaRF,...
				{SCgetOrds(SC.RING,'QF:'),'normal','individual',2E-3},... % {Ordinates, normal/skew, individual/family, deltaK}
				{SCgetOrds(SC.RING,'QD:'),'normal','individual',1E-3},...  % {Ordinates, normal/skew, individual/family, deltaK}
				{SCgetOrds(SC.RING,'QFA:'),'normal','family',4E-4},...  % {Ordinates, normal/skew, individual/family, deltaK}
				{SC.ORD.SkewQuad,'skew','individual',8E-3});  % {Ordinates, normal/skew, individual/family, deltaK}
		end	

		% Run LOCO
		[~,BPMData, CMData, FitParameters, LocoFlags, LOCOmodel] = loco(LocoMeasData,  BPMData,  CMData,  FitParameters,  LocoFlags,  LOCOmodel);
		
		% Apply lattice correction step
		SC = SClocoLib('applyLatticeCorrection',SC,FitParameters);
		
		% Apply orbit correction
		SC = SClocoLib('applyOrbitCorrection',SC);
		
		% Update progress structure
		CorrProgress{end+1} = SC.RING;
		
		% Include chromaticity correction after 3 iterations
		if nIter > 3
			% Match chromaticity
			SC = SClocoLib('fitChromaticity',SC,SCgetOrds(SC.RING,{'SF:','SD:'}),'targetChrom',[1.2 1.1],'verbose',1);
			
			% Update progress structure
			CorrProgress{end+1} = SC.RING;
		end
		
		% Plot progress
		SClocoLib('plotStatus',SC,Init,BPMData,CMData) ;	
	end
	
	% Apply diagnostic correction
	SC = SClocoLib('applyDiagnosticCorrection',SC,CMstep,CMData,BPMData);
	
	% Plot LOCO process
	plotProgress(SC,CorrProgress,Init)

end

% Plot correction progress (for debugging)
function plotProgress(SC,CorrProgress,Dist)
	
	nSteps = length(CorrProgress);
	nPlots = 5;
	tmp=[];
	for nStep=1:nSteps
		tmpSC = SC;
		tmpSC.RING = CorrProgress{nStep};
		OUT = SCcalcLatticeProperties(tmpSC);
		tmp.betaBeat(nStep,:)     = OUT.betaBeat;
		tmp.dispBeat(nStep,:)     = OUT.dispBeat;
		tmp.tuneShift(nStep,:)    = OUT.tuneShift;
		tmp.Chromaticity(nStep,:) = OUT.Chromaticity;
		tmp.Coupling(nStep,:)     = OUT.Coupling;
		tmp.Emittance(nStep,:)    = OUT.Emittance;
		tmp.orbit(nStep,:)        = OUT.orbit;
	end
	titleBBstr  = {'hor. \beta-beat','ver. \beta-function'};titleDBstr = {'hor. \Delta\eta','ver. \Delta\eta'};
	
	figure(446),clf
	for nDim=1:2
		subplot(2,nPlots,nPlots*(nDim-1)+1)
		plot(1E2*repmat(Dist.betaBeat(nDim),nSteps,1),'k','LineWidth',2);hold on
		plot(1E2*tmp.betaBeat(:,nDim),'O--','LineWidth',2);
		xlabel('Iteration');ylabel('\Delta\beta/\beta_0 [%]');title(titleBBstr{nDim})
		set(gca,'YScale','log')%,'YLim',yLim)
		
		subplot(2,nPlots,nPlots*(nDim-1)+2)
		plot(1E3*repmat(Dist.dispBeat(nDim),nSteps,1),'k','LineWidth',2);hold on
		plot(1E3*tmp.dispBeat(:,nDim),'O--','LineWidth',2);
		xlabel('Iteration');ylabel('\Delta\eta [mm]');title(titleDBstr{nDim})
		
		subplot(2,nPlots,nPlots*(nDim-1)+3)
		plot(abs(repmat(Dist.tuneShift(nDim),nSteps,1)),'k','LineWidth',2);hold on
		plot(abs(tmp.tuneShift(:,nDim)),'O--','LineWidth',2);
		xlabel('Iteration');ylabel('\Delta\nu [mm]');title('Tune')
		set(gca,'YScale','log')
		
		subplot(2,nPlots,nPlots*(nDim-1)+4)
		plot(repmat(Dist.Chromaticity(nDim),nSteps,1),'k','LineWidth',2);hold on
		plot(tmp.Chromaticity(:,nDim),'O--','LineWidth',2);
		xlabel('Iteration');ylabel('\xi');title('Chromaticity')
		
		subplot(2,nPlots,nPlots*(nDim-1)+5)
		plot(repmat(Dist.Coupling,1,nSteps)','k','LineWidth',2);hold on
		plot(tmp.Coupling,'O--','LineWidth',2);
		xlabel('Iteration');ylabel('RM');title('Coupling')
		set(gca,'YScale','log')
	end
	set(findall(gcf,'-property','FontSize'),'FontSize',20);set(gcf,'color','w');drawnow
	
end
