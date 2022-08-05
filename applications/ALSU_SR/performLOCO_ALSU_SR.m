function SC = performLOCO_ALSU_SR(SC,BPMords,CMords,varargin)
% performLOCO_ALSU_SR
% ===================
%
% NAME
% ----
% performLOCO_ALSU_SR - Performs LOCO for the ALSU-SR
%
% SYNOPSIS
% --------
% `SC = performLOCO_ALSU_SR(SC,BPMords,CMords [,options])`
%
% DESCRIPTION
% -----------
% This function performs LOCO for the ALSU-SR lattice.
%
% INPUT
% -----
% `SC`::
% 	The SC base structure
% `BPMords`::
% 	BPM ordinates used for orbit correction.
% `CMords`::
% 	CM ordinates used for orbit correction
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'nIter'` (1:8)`)::
%	Iterations to be carried out.
% `'CMstep'` (50E-6)::
% 	CM steps for RM measurement.
% `'deltaRF'` (1500)::
% 	RF step for dispersion measurement.
% `'RMpar'` ({})::
% 	Parameters for RM measurement.
% `'HorizontalDispersionWeight'` (1)::
%	Horizontal dispersion weight for LOCO.
% `'VerticalDispersionWeight'` (1)::
%	Vertical dispersion weight for LOCO.
% `'measRMeveryStep'` (0)::
%	Specify if RM should be re-measured after each correction step.
% `'RMstruct'` ([])::
%	Structure containing all information relevant for orbit correction.
%
% RETURN VALUE
% ------------
% `SC`::
% 	The SC base structure with corrected lattice.
%
% SEE ALSO
% --------
% *performDiagCorr_ALSU_SR*, *performOrbitCorr_ALSU_SR*
	
	
	% Parse optional arguments
	p = inputParser;
	addOptional(p,'nIter',1:8);
	addOptional(p,'RMpar',{});
	addOptional(p,'diagCorr',1);
	addOptional(p,'CMstep',50E-6); % CM step for response matrix calculation
	addOptional(p,'deltaRF',1500); % RF step for dispersion calculation
	addOptional(p,'RMfolder','~/sc/applications/ALSU_SR/LOCO/data/210106linear/'); % 
	addOptional(p,'RMstruct',[]); % RM used for orbit correction (inkl. dispersion)
	addOptional(p,'HorizontalDispersionWeight',1); % Hor. dispersion wight in LOCO (not used in diag-correction)
	addOptional(p,'VerticalDispersionWeight',1); % Ver. dispersion wight in LOCO (not used in diag-correction)
	parse(p,varargin{:});
	par = p.Results;
	
	% For use of precalculated LOCO RMs
	metaData.RMfolder   = par.RMfolder;
	metaData.usePreCalc = 0;
		
	% Load CM calibration offsets
	load('210103_SR_v21_CM_cal_Offsets');

	% Get orbit response matrix and dispersion from measurement
	LocoMeasData = SClocoLib('getMeasurement',SC,par.CMstep,par.deltaRF,SC.ORD.BPM,SC.ORD.CM,par.RMpar{:});

	% Set up BPM and CM data structure structures
	[BPMData,CMData] = SClocoLib('getBPMCMstructure',SC,par.CMstep);
 
	% Setup LOCO model
	[RINGdata,LocoFlags,Init] = SClocoLib('setupLOCOmodel',SC,...
		'HorizontalDispersionWeight',par.HorizontalDispersionWeight,...
		'VerticalDispersionWeight'  ,par.VerticalDispersionWeight,...
		'SVmethod',1E-2);
	
	% Define CM calibration offsets
	fields={'H','V'};
	for nDim=1:2
		i=1;
		for ord=CMData.([fields{nDim} 'CMIndex'])'
			FamInd = ~cellfun(@ isempty,regexp(SC.RING{ord}.FamName,CMfams{nDim}));
			CMcalOffsets{nDim}(i,1) = CMcalOffsetFam{nDim}(FamInd);
			i = i +1;
		end
		% Apply CM calibration offsets
		CMData.([fields{nDim} 'CMKicks']) = (1 + CMcalOffsets{nDim}) .* CMData.([fields{nDim} 'CMKicks']);
	end
	
	
	% Initialize correction progress structure (just for plotting)
	CorrProgress.RING{1} = SC.RING;
	CorrProgress.PAR{1}  = [];
		
	
	% Calculate model response matrix
	if isempty(par.RMstruct)
		par.RMstruct.RM         = SCgetModelRM(SC,BPMords,CMords,'trackMode','ORB','useIdealRing',1);
		par.RMstruct.eta        = SCgetModelDispersion(SC,BPMords,SC.ORD.Cavity);
		par.RMstruct.scaleDisp  = 1E7;
		par.RMstruct.BPMords    = BPMords;
		par.RMstruct.CMords     = CMords;
		par.RMstruct.alpha      = 5;
		par.RMstruct.CMsteps{1} = 50E-6 ./ mean(abs(RMstruct.RM(1:length(BPMords),1:length(CMords{1}))));
		par.RMstruct.CMsteps{2} = 50E-6 ./ mean(abs(RMstruct.RM(length(BPMords)+1:end,length(CMords{1})+1:end)));
	end
	
	
	%% Run LOCO
	for nIter=par.nIter
			
		if nIter==1
			% Set up LOCO fit parameter structure
			FitParameters = SClocoLib('setupFitparameters',SC,Init.SC.RING,RINGdata,par.deltaRF,...
				{SCgetOrds(SC.RING,'^QF1|^QD1'),'normal','individual',1E-2*13});   % {Ordinates, normal/skew, individual/family, deltaK}
		elseif nIter==2
			% Set up LOCO fit parameter structure
			FitParameters = SClocoLib('setupFitparameters',SC,Init.SC.RING,RINGdata,par.deltaRF,...
				{SCgetOrds(SC.RING,'^QF1|^QD1')     ,'normal','individual',1E-2*13},...
				{SCgetOrds(SC.RING,'^QF4|^QF5|^QF6'),'normal','individual',1E-2*20});   % {Ordinates, normal/skew, individual/family, deltaK}
					
		elseif nIter==3
			% Set up LOCO fit parameter structure
			FitParameters = SClocoLib('setupFitparameters',SC,Init.SC.RING,RINGdata,par.deltaRF,...
										{SCgetOrds(SC.RING,'^QF1|^QD1')     ,'normal','individual',1E-2*13},...
										{SCgetOrds(SC.RING,'^QF4|^QF5|^QF6'),'normal','individual',1E-2*20},...
										{SC.ORD.SkewQuad                    ,'skew'  ,'individual',100E-3});   % {Ordinates, normal/skew, individual/family, deltaK}
			
			LocoFlags.Coupling   = 'Yes'; % Include off-diagonal ORM elements
			LocoFlags.Dispersion = 'Yes'; % Include dispersion
		
		elseif nIter==5
		
			if par.diagCorr
				% Correct diagnostic devices
				SC = performDiagCorr_ALSU_SR(SC,SC.ORD.BPM,SC.ORD.CM,par.CMstep,par.deltaRF,par.RMpar,...
												'HorizontalDispersionWeight',par.HorizontalDispersionWeight,...
												'VerticalDispersionWeight'  ,par.VerticalDispersionWeight);
				
				% Setup LOCO model
				[RINGdata,LocoFlags,Init] = SClocoLib('setupLOCOmodel',SC,...
														'HorizontalDispersionWeight',par.HorizontalDispersionWeight,...
														'VerticalDispersionWeight'  ,par.VerticalDispersionWeight,...
														'SVmethod',1E-2);
				
				% Set up BPM and CM data structure structures
				[BPMData,CMData] = SClocoLib('getBPMCMstructure',SC,par.CMstep);
				% Include calibration offsets
				for nDim=1:2
					CMData.([fields{nDim} 'CMKicks']) = (1 + CMcalOffsets{nDim}) .* CMData.([fields{nDim} 'CMKicks']);
				end
				
				% Get orbit response matrix and dispersion from measurement
				LocoMeasData = SClocoLib('getMeasurement',SC,par.CMstep,par.deltaRF,SC.ORD.BPM,SC.ORD.CM,par.RMpar{:});
			end
		
			% Set up LOCO fit parameter structure
			FitParameters = SClocoLib('setupFitparameters',SC,Init.SC.RING,RINGdata,par.deltaRF,...
				{SCgetOrds(SC.RING,'^QF1|^QD1')     ,'normal','individual',1E-2*13},...
				{SCgetOrds(SC.RING,'^QF4|^QF5|^QF6'),'normal','individual',1E-2*20},...
				{SC.ORD.SkewQuad                    ,'skew'  ,'individual',70E-3});   % {Ordinates, normal/skew, individual/family, deltaK}
			
			LocoFlags.Coupling   = 'Yes'; % Include off-diagonal ORM elements
			LocoFlags.Dispersion = 'Yes'; % Include dispersion

		elseif nIter==7
			
			% Set up LOCO fit parameter structure
			FitParameters = SClocoLib('setupFitparameters',SC,Init.SC.RING,RINGdata,par.deltaRF,...
				{SCgetOrds(SC.RING,'^QF1|^QD1')     ,'normal','individual',1E-2*13},...
				{SCgetOrds(SC.RING,'^QF4|^QF5|^QF6'),'normal','individual',1E-2*20},...
				{SCgetOrds(SC.RING,'BEND2|BEND3')   ,'normal','individual',1E-2*7},...
				{SC.ORD.SkewQuad                    ,'skew'  ,'individual',70E-3});   % {Ordinates, normal/skew, individual/family, deltaK}
			
		end
		
		
		% Set the right flags for loading the precalculated RMs
		switch LocoFlags.Coupling
			case 'Yes'
				metaData.coupl = 1;
			case 'No'
				metaData.coupl = 0;
		end
		switch LocoFlags.Dispersion
			case 'Yes'
				metaData.disp = 1;
			case 'No'
				metaData.disp = 0;
		end
		
		% Run LOCO
		[~,BPMData, CMData, FitParameters, LocoFlags, RINGdata] = locoTH(LocoMeasData,  BPMData,  CMData,  FitParameters,  LocoFlags,  RINGdata, metaData);
		
		% Apply lattice correction step
		SC = SClocoLib('applyLatticeCorrection',SC,FitParameters);
		
 		% Apply orbit correction
		[SC,~] = performOrbitCorr_ALSU_SR(SC,par.RMstruct);
		
		% Update progress structure
		CorrProgress.RING{end+1} = SC.RING;
		CorrProgress.PAR{end+1} = FitParameters;
		
		% Plot progress
		SClocoLib('plotStatus',SC,Init,BPMData,CMData) ;
		
		if nIter > 3
			% Match chromaticity
			SC = SClocoLib('fitChromaticity',SC,SCgetOrds(SC.RING,{'SF','SD'}),'targetChrom',[2 1],'verbose',1);

			% Apply orbit correction
			[SC,~] = performOrbitCorr_ALSU_SR(SC,par.RMstruct);

			% Update progress structure
			CorrProgress.RING{end+1} = SC.RING;
			
		end
	end

% 	plotProgress(SC,CorrProgress,Init)
	
end
	
	
	

	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot correction progress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotProgress(SC,CorrProgress,Dist)
			
	OUT = calcLatticeProperties_ALSU_SR(Dist.SC);
	Dist.betaBeat = OUT.betaBeat;
	
	nSteps = length(CorrProgress.RING);
	nPlots = 5;
	tmp=[];
	for nStep=1:nSteps
		tmpSC = SC;
		tmpSC.RING = CorrProgress.RING{nStep};
		OUT = calcLatticeProperties_ALSU_SR(tmpSC);
		tmp.betaBeat(nStep,:)     = OUT.betaBeat;
		tmp.dispBeat(nStep,:)     = OUT.dispBeat;
		tmp.tuneShift(nStep,:)    = OUT.tuneShift;
		tmp.Chromaticity(nStep,:) = OUT.Chromaticity;
		tmp.Coupling(nStep,:)     = OUT.Coupling;
		tmp.Emittance(nStep,:)    = OUT.Emittance;
		tmp.orbit(nStep,:)        = OUT.orbit;
	end
	% Fit parameters
	par=[];ords=[];
	for nStep=2:length(CorrProgress.PAR)
		for nPar=1:length(CorrProgress.PAR{nStep}.Params)
			ord = CorrProgress.PAR{nStep}.Params{nPar}.ElemIndex;
			par(nStep,ord) = CorrProgress.PAR{nStep}.Values(nPar)/CorrProgress.PAR{nStep}.OrigValues(nPar);
		end
	end
	par(par==0)=nan;
	figure(4467)
	plot(par)
	
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
