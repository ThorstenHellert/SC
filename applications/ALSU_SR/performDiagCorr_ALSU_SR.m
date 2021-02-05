function SC = performDiagCorr_ALSU_SR(SC,BPMords,CMords,CMstep,deltaRF,RMpar,varargin)
% performDiagCorr_ALSU_SR
% =======================
%
% NAME
% ----
% performDiagCorr_ALSU_SR - Performs LOCO-based correction of BPM and CM calibration factors
%
% SYNOPSIS
% --------
% `SC = performDiagCorr_ALSU_SR(SC,BPMords,CMords,CMstep,deltaRF,RMpar [,options])`
%
% DESCRIPTION
% -----------
% This function performs LOCO-based correction of BPM and CM calibration factors for the ALSU-SR
% lattice.
%
% INPUT
% -----
% `SC`::
% 	The SC base structure
% `BPMords`::
% 	BPM ordinates.
% `CMords`::
% 	CM ordinates.
% `CMstep`::
% 	CM steps for RM measurement.
% `deltaRF`::
% 	RF step for dispersion measurement.
% `RMpar`::
% 	Parameters for RM measurement.
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'nIter'` (1:5)`)::
%	Iterations to be carried out.
% `'HorizontalDispersionWeight'` (1)::
%	Horizontal dispersion weight for LOCO.
% `'VerticalDispersionWeight'` (1)::
%	Vertical dispersion weight for LOCO.
% `'measRMeveryStep'` (0)::
%	Specify if RM should be re-measured after each correction step.
%
% RETURN VALUE
% ------------
% `SC`::
% 	The SC base structure with corrected BPM/CM calibration errors.
%
% SEE ALSO
% --------
% *performLOCO_ALSU_SR*, *performOrbitCorr_ALSU_SR*

	
	
	
	p = inputParser;
	addOptional(p,'nIter',1:5);
	addOptional(p,'HorizontalDispersionWeight',1); % Hor. dispersion wight in LOCO (not used in diag-correction)
	addOptional(p,'VerticalDispersionWeight',1); % Ver. dispersion wight in LOCO (not used in diag-correction)
	addOptional(p,'measRMeveryStep',0); % Measure response matrix at each step
	parse(p,varargin{:});
	par = p.Results;
		
	% Use precalculated LOCO RMs
	metaData.usePreCalc = 0;

	for nIter=par.nIter

		if par.measRMeveryStep || nIter==1 || length(par.nIter)==1
			% Get orbit response matrix and dispersion from measurement
			[LocoMeasData,CMsteps]= SClocoLib('getMeasurement',SC,CMstep,deltaRF,BPMords,CMords,RMpar{:});
		end
		
		% Setup LOCO model
		[LOCOmodel,LocoFlags,Init] = SClocoLib('setupLOCOmodel',SC,...
			'HorizontalDispersionWeight',par.HorizontalDispersionWeight,...
			'VerticalDispersionWeight',par.VerticalDispersionWeight,...
			'SVmethod',1E-2,...
			'ClosedOrbitType','FixedMomentum',...
			'ResponseMatrixCalculator','full');
		

		% Set up BPM and CM data structure structures
		[BPMData,CMData] = SClocoLib('getBPMCMstructure',SC,CMsteps,...
			{'BPMords',BPMords,BPMords},...
			{'CMords',CMords{1},CMords{2}});
				
		% Fit calibration factors
		if nIter==1 || nIter==3
			BPMData.FitGains    = 'No';
			CMData.FitKicks     = 'Yes';
			LocoFlags.SVmethod  = 1E-3;
		elseif nIter==2 || nIter==4
			BPMData.FitGains    = 'Yes';
			CMData.FitKicks     = 'No';
			LocoFlags.SVmethod  = 1E-3;
		else
			BPMData.FitGains    = 'Yes';
			CMData.FitKicks     = 'Yes';
			LocoFlags.SVmethod  = 1E-4;
		end
		
		
		
		% Set up LOCO fit parameter structure
		FitParameters = SClocoLib('setupFitparameters',SC,Init.SC.RING,LOCOmodel,deltaRF);
		
		% Run LOCO
		[~,BPMData, CMData, ~, ~, ~] = locoTH(LocoMeasData,  BPMData,  CMData,  FitParameters,  LocoFlags,  LOCOmodel , metaData);
		
		% Apply diagnostic correction
		SC = SClocoLib('applyDiagnosticCorrection',SC,CMsteps,CMData,BPMData,'meanToZero',0);
		

		% Adjust measured response matrix manually for CM and BPM calibration errors
		if ~par.measRMeveryStep
			if strcmp(CMData.FitKicks,'Yes')
				% Get fitted CM calibration factors
				fitCalCM = [CMData.HCMKicks./CMsteps{1};CMData.VCMKicks./CMsteps{2}]'/1000/2;
				
				% Remove from measured response matrix
				LocoMeasData.M = LocoMeasData.M ./ repmat(fitCalCM,size(LocoMeasData.M,1),1);
			end
			if strcmp(BPMData.FitGains,'Yes')
				% Get fitted BPM calibration factors
				fitCalBPM  = [BPMData.HBPMGain;BPMData.VBPMGain];
				
				% Remove from measured response matrix
				LocoMeasData.M = LocoMeasData.M ./ repmat(fitCalBPM,1,size(LocoMeasData.M,2));
			end
		end

		% Apply orbit correction
		if strcmp(CMData.FitKicks,'Yes') && par.measRMeveryStep
			SC = SClocoLib('applyOrbitCorrection',SC,'alpha',1,'BPMords',BPMords,'CMords',CMords);
		end
	end
end
	
