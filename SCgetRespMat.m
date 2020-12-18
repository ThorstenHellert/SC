function [RM,Err,CMsteps] = SCgetRespMat(SC,Amp,BPMords,CMords,varargin)
% SCgetRespMat
% ============
%
% NAME
% ----
% SCgetRespMat - Simulates the response matrix measurement
%
% SYNOPSIS
% --------
% `[M, Err, Cmsteps] = SCgetRespMat(SC, Amp, BPMords, CMords [, options])`
%
%
% DESCRIPTION
% -----------
% Gets the beam based response matrix based on the current injection pattern and using the BPMs and
% CMs specified in `BPMords` and `CMords`, respectively. By default the response matrix is measrued
% using a fixed initial kick amplitude `Amp` (see options below), which can be given either as a 
% single value for all CMs or as a cell array with amplitudes for each CM specified in `CMords`. 
%
% If the applied CM setpoint was clipped because of CM limits, the measurement is repeated with the
% different direction, thus `-Amp`. If the beam transmission for any applied CM step is less then
% for the reference measurement, the CM step is scaled to 90% and the measurument is repeated.
%
% Additionally, the measurement mode can be specified as `fixedOffset`. In this case `Amp` specifies
% the maximum BPM offset difference which should be achieved for each CM step. The CM step is
% iterated three times for every CM to reach the desired change of BPM readings.
%
% INPUTS
% ------
% `SC`::
%	SC base structure
% `Amp`::
%	Amplitude of response matrix measurement step [m or rad], either single value or cell array 
%   defining the amplitude for every CM specified in `CMords`
% `BPMords`::
%	List of BPM ordinates at which the reading should be returned
% `CMords`::
%	List of CM ordinates at which the reading should be returned
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'mode'` (`'fixedKick'`):: 
%    Measurement mode, either `'fixedKick'` or `'fixedOffset'`
% `'nSteps'` (2):: 
%    Number of CM steps (1st CM step is considered the reference). If more than 2 steps are 
%    specified, the measurement is bi-directional
% `'fit'` (`'linear'`):: 
%    Fit method, either `'linear'` or `'quadratic'`
% `'verbose'` (0)::
%	If true, debug information is printed.
%
% RETURN VALUES
% -------------
% `M`:: Response matrix [m/rad]
% `Err`:: Chi squared errror of RM entries [m/rad]
% `CMsteps`:: Maximum CM steps used for RM measurement [rad]
%
% SEE ALSO
% --------
% *SCgetModelRM*, *SCgetBPMreading*, *SCsetCMs2SetPoints*, *SCgetCMSetPoints*

	
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Initialization
	
	% Parse optional arguments
	p = inputParser;
	addOptional(p,'mode','fixedKick');
	addOptional(p,'nSteps',2);
	addOptional(p,'fit','linear');
	addOptional(p,'verbose',0);
	parse(p,varargin{:});
	par = p.Results;
	
	% Quadratic fit doesn't work with two points
	if par.nSteps==2 && ~strcmp(par.fit,'linear')
		warning('Only linear fit method reasonable with two points.')
		par.fit = 'linear';
	end
	
	% Check if CM step is specified correctly
	if (~iscell(Amp) && ~length(Amp)==1 ) || (iscell(Amp) && ( length(Amp{1})~=length(CMords{1}) && length(Amp{2})~=length(CMords{2}) ))
		error('RM amplitude must be defined as single value or cell array matching the number of used HCM and VCM.')
	end
	
	% Expand amplitude for further evaluation
	if ~iscell(Amp)
		Amp = {repmat(Amp,length(CMords{1}),1),repmat(Amp,length(CMords{2}),1)};
	end
	
	% Switch orbit or tracking mode
	if strcmp(SC.INJ.trackMode,'ORB')
		par.nTurns = 1;
		if par.verbose;fprintf('Calculate orbit response matrix for %d BPMs and %d|%d CMs with mode ''%s'' and amplitde %.0e|%.0e using %d steps ...',length(BPMords),length(CMords{1}),length(CMords{2}),par.mode,mean(Amp{1}),mean(Amp{2}),par.nSteps);end
	else
		par.nTurns = SC.INJ.nTurns;
		if par.verbose;fprintf('Calculate %d-turn trajectory response matrix for %d BPMs and %d|%d CMs with mode ''%s'' and amplitde %.0e|%.0e using %d steps ...',SC.INJ.nTurns,length(BPMords),length(CMords{1}),length(CMords{2}),par.mode,mean(Amp{1}),mean(Amp{2}),par.nSteps);end
	end
	
	% Alocate structures
	RM      = nan(2*par.nTurns*length(BPMords),length(CMords{1})+length(CMords{2}));
	Err     = nan(2*par.nTurns*length(BPMords),length(CMords{1})+length(CMords{2}));
	CMsteps = {zeros(par.nSteps,length(CMords{1})),zeros(par.nSteps,length(CMords{2}))};
	
	% Calculate reference trajectory
	Bref      = reshape(SCgetBPMreading(SC,'BPMords',BPMords)',[],1);
	
	if strcmp(SC.INJ.trackMode,'ORB') && any(isnan(Bref))
		error('No closed orbit found.')
	end
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Main function
	
	
	
	i = 1;
	for nDim=1:2
		
		% Get all used CM initial setpoints
		cmstart = SCgetCMSetPoints(SC,CMords{nDim},nDim);
		
		% Loop over CMs
		for nCM = 1:length(CMords{nDim})
			
			% Get maximum kick amplitude for CM based on method and beam reach
			[MaxStep,dB] = getKickAmplitude(SC,Bref,BPMords,CMords{nDim}(nCM),Amp{nDim}(nCM),nDim,par);
			
			% Get CM kick vector
			CMstepVec = linspace(-MaxStep,MaxStep,par.nSteps);
			
			% Check if only two CM settings are concernded : measurement already done by getKickAmplitude()
			if par.nSteps~=2
				
				% Initialize CM setpoint vector and BPM readings
				realCMsetPoint = cmstart(nCM) + CMstepVec;
				dB             = [zeros(length(CMstepVec)-1,length(Bref));dB'];

				% Loop over CM steps
				for nStep = 1:length(CMstepVec)
					
					% Don't evaluate reference point and point which was evaluated by getKickAmplitude() again
					if CMstepVec(nStep)~=0 && CMstepVec(nStep)~=MaxStep
						
						% Set corrector magnet
						[SC,realCMsetPoint(nStep)] = SCsetCMs2SetPoints(SC,CMords{nDim}(nCM),  cmstart(nCM) + CMstepVec(nStep)  ,nDim);
						
						% Calculate BPM reading differences
						dB(nStep,:) = reshape(SCgetBPMreading(SC,'BPMords',BPMords)',[],1) - Bref;
					end
				end
				
				% Get real CM differences
				dCM = realCMsetPoint-cmstart(nCM);
			else
				% Get CM differences
				dCM = MaxStep;
			end		
			
			CMsteps{nDim}(:,nCM) = dCM;
			
			% Calculate response matrix elements
			if par.nSteps==2
				RM(:,i) = dB / dCM;
			else
				for nBPM=1:size(dB,2)
					x = dCM(~isnan(dB(:,nBPM)))';
					y = dB(~isnan(dB(:,nBPM)),nBPM);
					switch par.fit
						case 'linear'
							RM(nBPM,i) = x\y;
						case 'quadratic'
							tmp = polyfit(x,y,2);
							RM(nBPM,i) = tmp(2);
					end
					% Calculate response matrix element error
					Err(nBPM,i) = sqrt(mean( (RM(nBPM,i)*dCM(~isnan(dB(:,nBPM)))-dB(~isnan(dB(:,nBPM)),nBPM)').^2 ));
				
				end
			end
			
			% Increase RM column
			i = i+1;
			
			% Reset corrector magnet
			[SC,~] = SCsetCMs2SetPoints(SC,CMords{nDim}(nCM),cmstart(nCM),nDim);
			
		end
	end
	
	
	% Set nan in response matrix to zero
	RM(isnan(RM)) = 0;
	
	if par.verbose;fprintf(' done.\n');end
	
end


% Get maximum kick amplitude for CM based on method and beam reach
function [MaxStep,dB] = getKickAmplitude(SC,Bref,BPMords,CMord,Amp,nDim,par)
	
	% Get CM initial setpoint
	cmstart = SCgetCMSetPoints(SC,CMord,nDim);
	
	% Guess initial CMstep
	MaxStep = Amp;
	
	switch par.mode
		case 'fixedKick'
			% Try full amplitude, if beam gets lost reduce amplitude
			for n = 1:20
				
				% Set corrector magnet
				[SC,realCMsetPoint] = SCsetCMs2SetPoints(SC,CMord,  cmstart + MaxStep  ,nDim);
				
				% If CM stauration, try different direction
				if realCMsetPoint ~= (cmstart + MaxStep)
					if par.verbose; fprintf('CM  clipped. Using different CM direction.\n'); end
					MaxStep = - MaxStep;
					[SC,~] = SCsetCMs2SetPoints(SC,CMord,  cmstart + MaxStep  ,nDim);
				end
				
				% Calculate BPM reading
				B = reshape(SCgetBPMreading(SC,'BPMords',BPMords)',[],1);
				
				% Calculate beam dump
				maxpos    = min([find(isnan(B   ),1,'first')-1 , par.nTurns*length(BPMords) ]);
				maxposRef = min([find(isnan(Bref),1,'first')-1 , par.nTurns*length(BPMords) ]);
				
				% Check if beam gets lost earlier then at reference trajectory
				if ~(maxpos < maxposRef)
					dB = B - Bref; 
					% Exit loop
					break
				else
					% CM step scaling factor will be reduced
					MaxStep = 0.9*MaxStep;
					if par.verbose; fprintf('Insufficient beam reach (%d/%d). CMstep reduced to %.1furad.\n',maxpos,maxposRef,1E6*MaxStep);end
				end
			end
			
		case 'fixedOffset'
			
			% Iterate CM settings
			for n=1:4
				
				% Set corrector magnet
				[SC,realCMsetPoint] = SCsetCMs2SetPoints(SC,CMord,  cmstart + MaxStep  ,nDim);
				
				% If CM stauration, try different direction
				if realCMsetPoint ~= (cmstart + MaxStep)
					if par.verbose; fprintf('CM  clipped. Using different CM direction.\n'); end
					MaxStep = - MaxStep;
					[SC,~] = SCsetCMs2SetPoints(SC,CMord,  cmstart + MaxStep  ,nDim);
				end
								
				% Calculate BPM reading
				B = reshape(SCgetBPMreading(SC,'BPMords',BPMords)',[],1);
				
				% Calculate beam dump
				maxpos    = min([find(isnan(B   ),1,'first')-1 , par.nTurns*length(BPMords) ]);
				maxposRef = min([find(isnan(Bref),1,'first')-1 , par.nTurns*length(BPMords) ]);
				
				% Check if beam gets lost earlier then at reference trajectory
				if (maxpos < maxposRef)
					MaxStep = 0.5*MaxStep;
					% CM step scaling factor will be reduced
					if par.verbose; fprintf('Insufficient beam reach (%d/%d). CMstep reduced to %.1furad.\n',maxpos,maxposRef,1E6*MaxStep); end
					continue
				end
				
				% Scale CM step
				MaxStep = MaxStep * Amp / max(abs(B - Bref));
			end
			
			% Calculate final BPM reading
			dB = reshape(SCgetBPMreading(SC,'BPMords',BPMords)',[],1) - Bref;
					
% 			fprintf('RMS BPM diff = %.3e\n',sqrt(mean(dB.^2)))
% 			figure(1);clf
% 			plot(1E6*dB)
% 			drawnow
	end
	
end
