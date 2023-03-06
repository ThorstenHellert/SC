function [SC,errorFlags] = SCBBA(SC,BPMords,magOrds,varargin)
% SCBBA
% =====
%
% NAME
% ----
% SCBBA - Performs model independent beam based alignment
%
% SYNOPSIS
% --------
% `[SC, errorFlags] = SCBBA(SC, BPMords, magOrds, [, options])`
%
%
% DESCRIPTION
% -----------
% Perform a model independend beam based alignment procedure using either two-turn
% trajectories or closed orbit bumps.
% In two-turn mode, for each BPM the injected beam trajectory is varied within a
% user defined range and the corresponding magnet is exercised on a user
% defined vector of setpoints. If the initial trajectory variation causes beam
% losses, the trajectory variation is reduced until all injected beam
% trajectories reach the BBA-BPM. If the final trajectory variation at the
% considered BBA-BPM is below a user defined threshold, a quadrupole may be
% exercised to change the phase advance between the injection point and the
% BPM (see option `'quadOrdPhaseAdvance'`). 
% In orbit mode an orbit bump is generated at each considerded BPM using orbit feedback 
% with weighting factors on a user defined window around the BBA BPM.
% Finally, for each injected beam trajectory or orbit bump step the offset variation due
% the change in magnet strength (see option `'magOrder'`) is recorded at the BPMs used 
% for measurement. A polynomial fit (see option `'fitOrder'`) is used to determine the
% center of the BBA-BPM with respect to the used magnet.
%
% INPUTS
% ------
%
% `SC`::
%	SC base structure
% `BPMords`::
%	[2 x n] array of BPM ordinates
% `magOrds`::
%	[2 x n] array of magnet ordinates for the corresponding BPMs in `BPMords`
%
% RETURN VALUES
% -------------
%
% `SC`::
%	SC base structure with updated BPM offsets
% `errorFlags`::
%	[2 x n] array of error flags for each BPM
%		
%  
% OPTIONS
% -------
% The following options can be specified as name-value pairs:
%
% `'mode'` (`SC.INJ.trackMode`)::
%    Orbit or two-turn trajectory mode (`'TBT'`)
% `'outlierRejectionAt'` (`Inf`)::
%    If the calculated BPM offset change is above the specified value, the measurement
%    is discarted and the BPM offset is not updated.
% `'fakeMeasForFailures'` (`1`)::
%    This option intends to mimic the operater's ability to identify errors in the 
%    measurement procedure and adjust the fine tuning parameters for individual BPMs.
%	 After performing the measurement routine, the rms value of the difference between
%    the BPM offsets and the magnet centers is calculated for both planes for all 
%    successful BPMs. If this flag is set to `1`, all BPM offsets at which the measurement
%    failed are artificially generated using a Gaussian distribution with two sigma cutoff
%    with the rms value as described above.
% `'dipCompensation'` (`1`)::
%	 Flag specifying if dipole compensation at combined function magnets should be used in
%    measurement. This works only if the considered magnet is equipped with a registered HCM.
%    If so, the HCM is used to compensate the bending angle change when changing the main 
%    quadrupole coil. If this option is switched off, the HCM is not used and the 'real' 
%    quadrupole center is determined in the measurement and the BPM offset is adjusted in 
%    order to account for the difference between the real quadrupole center and its design 
%    value.
% `'nSteps'` (`10`)::
%	 Number of different trajectories/orbit bump steps for each magnet setting
% `'fitOrder'` (`1`)::
%	 Order of polynominal fitting, e.g. `1` for linear, `2` for quadratic fit.
% `'magOrder'` (`2`)::
%	 Order of magnet setpoint change, e.g. `2` for quadrupole, `3` for sextupole etc.
% `'magSPvec'` (`[0.95 1.05]`)::
%	 Strength variation of magnets. Can be either 1xN array of setpoint variations
%    (applied to all specified magnets) or cell array of 1xN arrays with
%    size(magSPvec)==size(magOrds) such that each magnet has it's individual setpoint
%    variation array.
% `'magSPflag'` (`'rel'`)::
%	 Specify how magnet setpoint as specified by option 'magSPvec' should be changed,
%    e.g. relative or absolute (see *SCsetMags2SetPoints*).
% `'skewQuadrupole'` (`0`)::
%	 If true, it is assumed that a skew quadrupole is used. Thus, the BPM readings in the 
%    dimension other than the trajectory/orbit excitation is used for evaluation.
% `'switchOffSext'` (`0`)::
%	 Flag specifying if sextupole coil in BBA magnet should be switched off 
%    (e.g. if quadrupole trim coils are used).
% `'RMstruct'` (`[]`)::
%	 Structure containing pre-calculated response matrix etc. for orbit feedback
%    (orbit mode only). If empty the relevant arrays are calculated with default
%    parameters.
% `'orbBumpWindow'` (`5`)::
%	 Number of BPMs adjacent to BBA BPM (upstream and downstream) where the BPM weighting 
%    in orbit feedback is set to zero to allow for a pseudo orbit bump (orbit mode only).
% `'useBPMreadingsForOrbBumpRef'` (`0`)::
%	 If true the actual BPM readings will be used for the reference when generating the
%    pseudo orbit bump (orbit mode only) instead of zeros. (See github issue #22)
% `'BBABPMtarget'` (`1E-3`)::
%	 Desired offset variation at BBA-BPM (BPM adjacent to magnet) which should be
%    achieved by trajectory change or orbit bump.
% `'minBPMrangeAtBBABBPM'` (`0.5E-3`)::
%	 Threshold of change of offset at BBA-BPM; if below BBA evaluation is not performed.
% `'minBPMrangeDownstream'` (`100E-6`)::
%	 Minimum change of offset at downstream BPMs to be included in calculation.
% `'nXPointsNeededAtMeasBPM'` (`3`)::
%	 Number of x-positions at BBA-BPM required for linear regression at
%	 measurement BPMs.
% `'maxNumOfDownstreamBPMs'` (`length(SC.ORD.BPM)`)::
%	 Number downstream BPMs considered in the data evaluation (2-turn mode only).
% `'minSlopeForFit'` (`0.03`)::
%	 Minimum fitted slope at measurement BPMs with respect to magnet change which
%    is still used in determining BPM center (a very small slope usually indicates
%    an unfortunate phase advance and can significantly affect the measurement if
%    not excluded). Only with linear fitting.
% `'maxStdForFittedCenters'` (`600E-6`)::
%	 If standard deviation of the fitted BPM centers as determined by all
%	 downstream BPMs exceeds this value, output will be `'nan'`.
% `'maxTrajChangeAtInjection'` (`[0.9E-3 0.9E-3]`)::
%	 Maximum offset [m] and kick [rad] change at injection (2-turn mode only).
% `'quadOrdPhaseAdvance'` (`[]`)::
%	 Magnet ordinate which is used to change the phase advance between injection
%    and BBA-BPM (2-turn mode only).
% `'quadStrengthPhaseAdvance'` (`[0.95 1.05]`)::
%	 Relative magnet strength variation used to change the phase advance between
%    injection and BBA-BPM (2-turn mode only).
% `'plotLines'` (0)::
%	 If true, each injected beam and intermediate BBA results will be plotted.
% `'plotResults'` (0)::
%    If true, final BBA results are plotted.
% `'verbose'` (0)::
%	 If true, debug information is printed.
%
% ERROR FLAGS
% -----------
% The [2 x n] array of error flags specify if and why the measurement failed for each BPM 
% and may have the following value:
%
% (0):: All good
% (1):: Max. range at BBA-BPM to small (see option 'minBPMrangeAtBBABBPM')
% (2):: Max. range at downstream BPM to small (see option 'minBPMrangeOtherBPM')
% (3):: Fitted magnetic centers to far spread out (see option 'maxStdForFittedCenters')
% (4):: All downstream BPM measurements failed
% (5):: Unexpected error during data evaluation
% (6):: Calculated BPM offset change too large (see option 'outlierRejectionAt')
%
% SEE ALSO
% --------
% *SCsetMags2SetPoints*, *SCgetBPMreading*, *SCfeedbackRun*

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Input check
	
	p = inputParser;
	addOptional(p,'mode',SC.INJ.trackMode);
	addOptional(p,'outlierRejectionAt',Inf);
	addOptional(p,'fakeMeasForFailures',0);
	addOptional(p,'dipCompensation',1);
	addOptional(p,'nSteps',10);
	addOptional(p,'fitOrder',1);
	addOptional(p,'magOrder',2);
	addOptional(p,'magSPvec',[0.95,1.05]);
	addOptional(p,'magSPflag','rel');
	addOptional(p,'skewQuadrupole',0);
	addOptional(p,'switchOffSext',0);
	addOptional(p,'RMstruct',[]);
	addOptional(p,'orbBumpWindow',5);
	addOptional(p,'useBPMreadingsForOrbBumpRef',0);
	addOptional(p,'BBABPMtarget',1E-3);
	addOptional(p,'minBPMrangeAtBBABBPM',500E-6);
	addOptional(p,'minBPMrangeOtherBPM',100E-6);
	addOptional(p,'maxStdForFittedCenters',600E-6);
	addOptional(p,'nXPointsNeededAtMeasBPM',3);
	addOptional(p,'maxNumOfDownstreamBPMs',length(SC.ORD.BPM));
	addOptional(p,'minSlopeForFit',0.03);
	addOptional(p,'maxTrajChangeAtInjection',[.9E-3 .9E-3]);
	addOptional(p,'quadOrdPhaseAdvance',[ ]);
	addOptional(p,'quadStrengthPhaseAdvance',[0.95 1.05]);
	addOptional(p,'plotLines',0);
	addOptional(p,'plotResults',0);
	addOptional(p,'verbose',0);
	parse(p,varargin{:});
	par = p.Results;

	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Initialization
	
	if any(size(BPMords)~=size(magOrds))
		error('Input arrays for BPMs and magnets must be same size.')
	end
	
	if strcmp(par.mode,'TBT') && SC.INJ.nTurns~=2
		fprintf('Setting number of turns to 2.\n')
		SC.INJ.nTurns = 2;
	end
	
	% Expand setpoint variation array
	if ~iscell(par.magSPvec)
		par.magSPvec = repmat({par.magSPvec},size(magOrds));
	end
	
	% Save initial BPM offset errors (for plotting)
	initOffsetErrors = getBPMoffsetFromMag(SC,BPMords,magOrds);
	
	% Initialize error flags
	errorFlags = nan(size(BPMords));
	
	% Offset and kick angle variation at injection (for trajectory mode)
	kickVec0  = par.maxTrajChangeAtInjection' .* repmat(linspace(-1,1,par.nSteps),2,1);

	% Save initial injection trajectory (for trajectory mode)
	initialZ0 = SC.INJ.Z0;
	
	% Calculate all relevant arrays for orbit correction if not supplied by user
	if strcmp(par.mode,'ORB')
		if isempty(par.RMstruct)
			fprintf('Calculating orbit response matrix and dispersion.\n')
			par.RMstruct.RM        = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'trackMode','ORB','useIdealRing',1); 
			par.RMstruct.eta       = SCgetModelDispersion(SC,SC.ORD.BPM,SC.ORD.Cavity);
			par.RMstruct.scaleDisp = 1E7;
			par.RMstruct.BPMords   = SC.ORD.BPM;
			par.RMstruct.CMords    = SC.ORD.CM;
			par.RMstruct.MinvCO    = SCgetPinv([par.RMstruct.RM par.RMstruct.scaleDisp*par.RMstruct.eta],'alpha',5);
		end
	end
	
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Perform BBA routine
	
	% Loop over BPMs
	for jBPM=1:size(BPMords,2) % jBPM: Index of BPM adjacent to magnet for BBA

		% Horizontal/vertical
		for nDim=1:size(BPMords,1)
			if par.verbose;fprintf('BBA-BPM %d/%d, nDim = %d\n',jBPM,size(BPMords,2),nDim);end

			% Save initial machine state (for convenience)
			SC0 = SC;
			
			% Get BPM index-ordinate pairing
			BPMind = find(BPMords(nDim,jBPM)==SC.ORD.BPM);
			
			% Define ordinate of magnet for BBA
			mOrd = magOrds(nDim,jBPM);		
			
			% Switch off sextupole coil at BBA magnet?
			if par.switchOffSext
				% Switch off sextupole coil
				SC = SCsetMags2SetPoints(SC,mOrd,2,3,0,'method','abs');
				% Run orbit feedback
				[SC,~] = SCfeedbackRun(SC,par.RMstruct.MinvCO,...
										'target',0,...
										'maxsteps',50,...
										'scaleDisp',par.RMstruct.scaleDisp,...
										'verbose',par.verbose,...
										'BPMords',par.RMstruct.BPMords,...
										'CMords',par.RMstruct.CMords,...
										'eps',1E-6);
			end

			
			switch par.mode
				case 'ORB'
					% Get orbit bump at BBA BPM
					[CMords,CMvec] = getOrbitBump(SC,mOrd,BPMords(nDim,jBPM),nDim,par);
					
					% Perform data measurement
					[BPMpos,tmpTra] = dataMeasurement(SC,mOrd,BPMind,jBPM,nDim,par,CMords,CMvec);

				case 'TBT'
					% Scale maximum initial trajectory variation before beam gets lost
					[kickVec, BPMrange] = scaleInjectionToReachBPM(SC,BPMind,nDim,initialZ0,kickVec0,par);
					
					% Check if offset range at BBA-BPM is to small and vary phase advance
					if ~isempty(par.quadOrdPhaseAdvance) && BPMrange < par.BBABPMtarget
						[SC,kickVec] = scanPhaseAdvance(SC,BPMind,nDim,initialZ0,kickVec0,par);
					end
					
					% Perform data measurement
					[BPMpos,tmpTra] = dataMeasurement(SC,mOrd,BPMind,jBPM,nDim,par,initialZ0,kickVec);
			end
			
			% Perform data evaluation
			[OffsetChange,errorFlags(nDim,jBPM)] = dataEvaluation(SC,BPMords,jBPM,BPMpos,tmpTra,nDim,mOrd,par);

			% Reset machine (should ideally be done by setting all magnet setpoints to intial values)
			SC = SC0;
						
			% Outlier rejection
			if  OffsetChange > par.outlierRejectionAt
				OffsetChange          = nan;
				errorFlags(nDim,jBPM) = 6;
			end
			
			% Check if BPM center could be identified in measurement
			if ~isnan(OffsetChange)
				% Apply BBA results
				SC.RING{BPMords(nDim,jBPM)}.Offset(nDim) = SC.RING{BPMords(nDim,jBPM)}.Offset(nDim) + OffsetChange;
			end
				
		end

		% Plot BBA results
		if par.plotResults
			plotBBAResults(SC,initOffsetErrors,errorFlags,jBPM,BPMords,magOrds)
		end
	end

	% Fake measurement to deal with failed BPMs
	if par.fakeMeasForFailures
		SC = fakeMeasurement(SC,BPMords,magOrds,errorFlags);
	end
% 	% Remove outliers
% 	SC = removeOutliers(SC,BPMords,magOrds,par);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data measurement
function [BPMpos,tmpTra] = dataMeasurement(SC,mOrd,BPMind,jBPM,nDim,par,varargin)
	% This function performs the data measurement by looping through the magnet
	% setpoints and changing either the injected beam trajectory or the CMs to create 
	% a local orbit bump.
	%
	% INPUTS
	% ------
	% `SC`:: SC base structure
	% `mOrd`:: 	ordinate of BBA magnet
	% `BPMind`:: Index of BBA BPM within SC.ORD.BPM
	% `nDim`:: Hor/Ver plane
	% `par`:: Auxiliary structure
	% `varargin{1}`:: Either CM ordinates or initial injected trajectory
	% `varargin{2}`:: Either CM setpoint matrix or injected trajectory offset and angle variations
	%
	% RETURN VALUES
	% -------------
	% `BPMpos`:: Offset variation at BBA BPM
	% `tmpTra`:: All other BPM readings throught the measurement

	% Check if skew quadrupole is used
    if par.skewQuadrupole
        magType = 1;
        if nDim==1
            measDim = 2;
		else
            measDim = 1;
        end
    else
        magType = 2;
        measDim = nDim;
    end
	
	
	switch par.mode
		case 'ORB'
			
			% Used CM ordinates and setpoint matrix
			CMords = varargin{1};
			CMvec  = varargin{2};
			
			% Measurement steps per magnet setting
			nMsteps = size(CMvec{nDim},1);

			% Prealocate x and y data
			tmpTra = nan(nMsteps,length(par.magSPvec{nDim,jBPM}),length(SC.ORD.BPM));
			BPMpos = nan(nMsteps,length(par.magSPvec{nDim,jBPM}));

		case 'TBT'
			% Initial trajectory and kick variations
			initialZ0 = varargin{1};
			kickVec   = varargin{2};
			
			% Measurment steps per magnet setting
			nMsteps = size(kickVec,2);
			
			% Prealocate x and y data
			tmpTra = nan(nMsteps,length(par.magSPvec{nDim,jBPM}),par.maxNumOfDownstreamBPMs);
			BPMpos = nan(nMsteps,length(par.magSPvec{nDim,jBPM}));
	end
	
	% Loop over different magnet settings
	for nQ=1:length(par.magSPvec{nDim,jBPM})

		% Set magnet to different setpoints
		SC = SCsetMags2SetPoints(SC,mOrd,magType,par.magOrder,par.magSPvec{nDim,jBPM}(nQ),...
			'method',par.magSPflag,...
			'dipCompensation',par.dipCompensation);

		% Loop over different trajectories
		for nKick=1:nMsteps

			switch par.mode
				case 'ORB'
					% Set CM
					for nD=1:2
						SC = SCsetCMs2SetPoints(SC,CMords{nD},CMvec{nD}(nKick,:),nD,'abs');
					end
				case 'TBT'
					% Vary initial trajectory
					SC.INJ.Z0(2*nDim)   = initialZ0(2*nDim  ) + kickVec(2,nKick); % kick angle
					SC.INJ.Z0(2*nDim-1) = initialZ0(2*nDim-1) + kickVec(1,nKick); % offset
			end
			
			% Calculate beam reading
			B = SCgetBPMreading(SC);

			% JUST FOR PLOTTING
			if par.plotLines; plotBBAstep(SC,BPMind,jBPM,nDim,nQ,mOrd,nKick,par);end

			% Save BBA-BPM readings
			BPMpos(nKick,nQ) = B(nDim,BPMind);
	
			switch par.mode
				case 'ORB'
					% Save all BPM readings
					tmpTra(nKick,nQ,:) = B(measDim, : );
				case 'TBT'
					% Save downstream BPM readings
					tmpTra(nKick,nQ,:) = B(measDim, (BPMind+1):(BPMind+par.maxNumOfDownstreamBPMs) );
			end
		end
	end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data evaluation
function [OffsetChange,Error] = dataEvaluation(SC,BPMords,jBPM,BPMpos,tmpTra,nDim,mOrd,par)
	% This function calculates the BPM offset from the measured data by identifying the zero 
	% crossing of each offset change. If 'par.plotLines=1', measurement details are plotted.
	% Note: it is recommended to comment out the line 'SCplotBPMreading(SC,B,T1);' in the 
	% function SCgetBPMreading() for performance reasons if detailed plotting is requested.
	%
	% INPUTS
	% ------
	% `SC`:: SC base structure
	% `BPMind`:: Index of BBA BPM within SC.ORD.BPM
	% `nDim`:: Hor/Ver plane
	% `initialZ0`:: Initial injected beam trajectory
	% `kickVec0`:: Initial injected beam trajectory change
	% `par`:: Auxiliary structure
	%
	% RETURN VALUES
	% -------------
	% `OffsetChange`:: Calculated BPM offset change (to be added on lattice element)
	% `Error`:: Flag specifying if and why measurement failed. Values are:
	%           (0): All good
	%           (1): Max. range at BBA-BPM to small (see option 'minBPMrangeAtBBABBPM')
	%           (2): Max. range at downstream BPM to small (see option 'minBPMrangeOtherBPM')
	%           (3): Fitted magnetic centers to far spread out (see option 'maxStdForFittedCenters')
	%           (4): All downstream BPM measurements failed 
	%           (5): Unexpected error
	
	% Plot measurement details
	if par.plotLines;figure(56),clf;hold on; p(1)=plot3(0,1E6*SC.RING{mOrd}.T2(2*nDim-1),0,'rD','MarkerSize',40,'MarkerFaceColor','b'); hold on;set(gca,'box','on');end
 	
	% Initialize output
	OffsetChange = nan;
	Error        = 5;
	
	% Initialize arrays
	tmpCenter = nan(1,(size(tmpTra,2)-1)*par.maxNumOfDownstreamBPMs);
	tmpNorm   = nan(1,(size(tmpTra,2)-1)*par.maxNumOfDownstreamBPMs);
	tmpRangeX = zeros(1,(size(tmpTra,2)-1)*par.maxNumOfDownstreamBPMs);
	tmpRangeY = zeros(1,(size(tmpTra,2)-1)*par.maxNumOfDownstreamBPMs);

	i = 0;
	% Loop over downstream BPMs
	for nBPM=1:par.maxNumOfDownstreamBPMs

		% y-data: downstream trajectory differences with respect to magnet setting
		y0 = diff(tmpTra(:,:,nBPM),1,2);
		% x data: BBA-BPM readings
		x0 = repmat(mean(BPMpos,2),1,size(y0,2));

		% Loop over different magnet differences
		for nKick=1:size(y0,2)
			i = i+1;

			% Current x and y data
			y = y0(:,nKick);
			x = x0(:,nKick);

			% Remove nan readings
			x(isnan(y)) = [];
			y(isnan(y)) = [];

			% Skip if BPM signal does not exist
			if isempty(x) || isempty(y)
				continue
			end

			% Calculate maximum trajectory change at BPMs
			tmpRangeX(i) = abs(min(x)-max(x));
			tmpRangeY(i) = abs(min(y)-max(y));

			% Initialize linear regression solution
			sol = nan(1,2);
			
			% Check if sufficient transmission and signal is achieved
			if length(x)>=par.nXPointsNeededAtMeasBPM ...           % Enough different points at BBA-BPM
					&& tmpRangeX(i) > par.minBPMrangeAtBBABBPM ...  % Range at BBA-BPM        > par.minBPMrangeAtBBABBPM
					&& tmpRangeY(i) > par.minBPMrangeOtherBPM       % Range at downstream BPM > pra.minBPMrangeOtherBPM

				if par.fitOrder==1
					% Linear regression
					sol = [ones(size(x)),x]\y;
					% Convert to polyfit convention (for plotting)
					sol = sol([2 1]);
					
					% Check if fitted slope is to small
					if abs(sol(1)) < par.minSlopeForFit
						sol(1) = nan;
					end
					
					% Calculate fitted magnetic center
					tmpCenter(i) = -sol(2)/sol(1);
					tmpNorm(i)   = 1/sqrt(sum((sol(1)*x+sol(2)-y).^2));
					
				else
					% Polynom fit
					[sol,S] = polyfit(x,y,par.fitOrder);
					
					% Calculate fitted magnetic center
					if par.fitOrder==2
						tmpCenter(i) = - (sol(2)/(2*sol(1)));
					else
						% Assume smallest root is the one
						tmpCenter(i) = min(abs(roots(sol)));
					end
					tmpNorm(i)   = 1/S.normr;
					

				end
				
% 				% Plot individual fit (helps with debugging)
% 				figure(42);clf;hold on
% 				plot(1E6*x,1E3*y,'x')
% 				plot(1E6*x,1E3*polyval(sol,x))
% 				plot(1E6*tmpCenter(i),1E3*polyval(sol,tmpCenter(i)),'kO')
% 				title(sprintf('Goodness of fit: %.1e',tmpNorm(i)))
% 				xlabel('BBA BPM offset [um]');ylabel('Downstream BPM offset change [mm]')
% 				legend({'Data','Fit','Root'});drawnow

			end

			
			% Just plotting
			if par.plotLines
				p(2)=plot3(repmat(jBPM+nBPM,size(x)),1E6*x,1E3*y,'ko');
				tmp=plot3(repmat(jBPM+nBPM,size(x)),1E6*x,1E3*polyval(sol,x),'k-');p(3)=tmp(1);
				p(4)=plot3(jBPM+nBPM,1E6*tmpCenter(nBPM),0,'Or','MarkerSize',10);
			end

		end
	end

	
	% Check if measurement requirements have been met
	if (max(tmpRangeX) < par.minBPMrangeAtBBABBPM)  
		% Max. range at BBA-BPM to small
		Error = 1;
	elseif max(tmpRangeY) < par.minBPMrangeOtherBPM 
		% Max. range at downstream BPM to small
		Error = 2;
	elseif std(tmpCenter,'omitnan') > par.maxStdForFittedCenters  
		% Fitted magnetic centers to far spread out
		Error = 3;
	elseif isempty(find(~isnan(tmpCenter),1))  
		% All downstream BPM measurements failed 
		Error = 4;
	else
		% Calculate mean magnetic center
		OffsetChange = sum(tmpCenter.*tmpNorm,'omitnan')/sum(tmpNorm,'omitnan');
		Error = 0;
	end

	% Check if dipole compensation is not used (for combined function quadrupoles)
	if ~par.dipCompensation && nDim==1 && SC.RING{mOrd}.NomPolynomB(2)~=0
		% Read magnet properties
		if isfield(SC.RING{mOrd},'BendingAngle')
			B = SC.RING{mOrd}.BendingAngle;
		else
			B = 0;
		end
		K = SC.RING{mOrd}.NomPolynomB(2);
		L = SC.RING{mOrd}.Length;
		
		% Adjust for design quadrupole offset
		OffsetChange = OffsetChange + B/L/K;
	end
	
	
	% Just plotting
	if par.plotLines
		p(5)=plot3(0,1E6*OffsetChange,0,'kD','MarkerSize',30,'MarkerFaceColor','r'); 
		p(6)=plot3(0,1E6*(SC.RING{BPMords(nDim,jBPM)}.Offset(nDim)+SC.RING{BPMords(nDim,jBPM)}.SupportOffset(nDim)+OffsetChange),0,'kD','MarkerSize',30,'MarkerFaceColor','g'); 
		title(sprintf('BBA-BPM: %d // mOrd: %d // mFam: %s // nDim: %d // FinOffset = %.0f $\\mu m$',jBPM,mOrd,SC.RING{mOrd}.FamName,nDim,1E6*abs(SC.RING{BPMords(nDim,jBPM)}.Offset(nDim) + SC.RING{BPMords(nDim,jBPM)}.SupportOffset(nDim) + OffsetChange - SC.RING{mOrd}.MagnetOffset(nDim) - SC.RING{mOrd}.SupportOffset(nDim))));
		legend(p,{'Magnet center','Measured offset change','Line fit','Fitted BPM offset (individual)','Fitted BPM offset (mean)','Predicted magnet center'})
		xlabel('Index of BPM'); ylabel('BBA-BPM offset [$\mu$m]'); zlabel('Offset change [mm]');
		set(findall(figure(56),'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');set(findall(gcf,'-property','FontSize'),'FontSize',18);set(gcf,'color','w');drawnow;
	end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale injection trajectory variation
function [kickVec, BPMrange] = scaleInjectionToReachBPM(SC,BPMind,nDim,initialZ0,kickVec0,par)
	% This function reduces the injected beam trajectory variation if the beam gets lost before
	% reaching the BBA BPM.	
	%
	% INPUTS
	% ------
	% `SC`:: SC base structure
	% `BPMind`:: Index of BBA BPM within SC.ORD.BPM
	% `nDim`:: Hor/Ver plane
	% `initialZ0`:: Initial injected beam trajectory
	% `kickVec0`:: Initial injected beam trajectory change
	% `par`:: Auxiliary structure
	%
	% RETURN VALUES
	% -------------
	% `kickVec`:: Potentially scaled injected beam trajectory change 
	% `BPMrange`:: Final offset range at BBA BPM

	% Initialize scaling factor and transmission falg
	fullTrans     = 0; % Transmission flag
	scalingFactor = 1; % Scaling factor for initial trajectory variation
	BPMrange      = 0; % Offset range at considered BPM

	% Loop as long as BPMs show nan readings and scaling factor is larger than numerical noise (zero)
	while fullTrans == 0 && scalingFactor > 1E-6
		% Reset lunch condition
		SC.INJ.Z0 = initialZ0;
		tmpBPMpos = [];

		% Check beam dump position when varying initial trajectory % % % % % % % % % % % % % % %
		for nK=[1 size(kickVec0,2)]
			% Vary initial trajectory
			SC.INJ.Z0(2*nDim)   = initialZ0(2*nDim  ) + scalingFactor * kickVec0(2,nK); % kick angle
			SC.INJ.Z0(2*nDim-1) = initialZ0(2*nDim-1) + scalingFactor * kickVec0(1,nK); % offset

			% Calculate beam reading
			B = SCgetBPMreading(SC);

			% Save relevant BPM reading
			tmpBPMpos(end+1) = B(nDim,BPMind);
		end

		% Check if BBA-BPM has NaN-reading
		if any(isnan(tmpBPMpos))
			% Reduce scaling factor
			scalingFactor = scalingFactor - 0.1;
		else
			% Set transmission flag to 1
			fullTrans = 1;
			BPMrange  = abs(max(tmpBPMpos)-min(tmpBPMpos));
		end
	end

	if scalingFactor<1E-6
		scalingFactor = 1;
		warning('Something went wrong. No beam transmission at all(?)')
	end

	kickVec       = scalingFactor .* kickVec0;

	if par.verbose;fprintf('Initial trajectory variation scaled to [%.0f|%.0f]%% of its initial value, BBA-BPM range %.0fum.\n',100*(kickVec([1 end])./kickVec0([1 end])),1E6*BPMrange);end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change BBA-BPM offset variation by a phase advance variation
function [SC,kickVec] = scanPhaseAdvance(SC,BPMind,nDim,initialZ0,kickVec0,par)
	% This function changes the phase advance by changing the strength of a quadrupole
	% in order to increase the offset variation at the BBA BPM (currently only in 2-turn mode).
	%
	% INPUTS
	% ------
	% `SC`:: SC base structure
	% `BPMind`:: Index of BBA BPM within SC.ORD.BPM
	% `nDim`:: Hor/Ver plane
	% `initialZ0`:: Initial injected beam trajectory
	% `kickVec0`:: Initial injected beam trajectory change
	% `par`:: Auxiliary structure
	%
	% RETURN VALUES
	% -------------
	% `SC`::  SC base structure with adjusted phase advance in lattice
	% `kickVec`:: Potentially scaled injected beam trajectory change 

	
	% Define quadrupole ordinate, strength vector and initial setpoint
	mOrd = par.quadOrdPhaseAdvance;
	qVec = par.quadStrengthPhaseAdvance;
	q0   = SC.RING{mOrd}.SetPointB(2);

	% Initialize all BPM ranges
	allBPMRange = zeros(length(qVec));

	% Loop over quadrupole setpoints
	for nQ=1:length(qVec)
		if par.verbose;fprintf('BBA-BPM range to small, try to change phase advance with quad ord %d to %.2f of nom. SP.\n',par.quadOrdPhaseAdvance,qVec(nQ));end

		% Set quadrupole to different setpoints
		SC = SCsetMags2SetPoints(SC,mOrd,2,2,qVec(nQ),'method','rel','dipCompensation',1);

		% Rescale initial trajectory variation
		[kickVec, BPMrange] = scaleInjectionToReachBPM(SC,BPMind,nDim,initialZ0,kickVec0,par);

		% Record BPM range
		allBPMRange(nQ) = BPMrange;

		if par.verbose;fprintf('Initial trajectory variation scaled to [%.0f|%.0f]%% of its initial value, BBA-BPM range %.0fum.\n',100*(kickVec([1 end])./kickVec0([1 end])),1E6*BPMrange);end

		% End loop if offset range at BBA-BPM is sufficient
		if ~( BPMrange < par.BBABPMtarget )
			break
		end
	end

	% Check if offset range at BBA-BPM is still to small
	if BPMrange < par.BBABPMtarget
		% Check if any BPM range is better than the initial one
		if BPMrange<max(allBPMRange)
			[~,nBest] = max(allBPMRange);
			% Set quadrupole to best setpoints
			SC = SCsetMags2SetPoints(SC,mOrd,2,2,qVec(nBest),'method','rel','dipCompensation',1);
			if par.verbose;fprintf('Changing phase advance of quad with ord %d NOT succesfull, returning to best value with BBA-BPM range = %.0fum.\n',mOrd,1E6*max(allBPMRange));end
		else
			% Reset quadrupole to initial setpoints
			SC = SCsetMags2SetPoints(SC,mOrd,2,2,q0,'method','abs','dipCompensation',1);
			if par.verbose;fprintf('Changing phase advance of quad with ord %d NOT succesfull, returning to initial setpoint.\n',mOrd);end
		end
	else
		if par.verbose;fprintf('Change phase advance of quad with ord %d successful. BBA-BPM range = %.0fum.\n',mOrd,1E6*BPMrange);end
	end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get orbit bump
function [CMords,CMvec] = getOrbitBump(SC,mOrd,BPMord,nDim,par)
	% This function is used to create a pseudo orbit bump by running orbit feedback on the actual 
	% machine with target 
	% of zero except at the BBA-BPM (here the target is 'par.BBABPMtarget'). The weights 
	% at the BPMs upstream and downstream the BBA BPM as defined by 'par.orbBumpWindow' are set to 
	% zero in order to give some slack.
	% Note that this could also be done by creating a 'real' orbit bump using 3 CMs, however we
	% found that this procedure is faster and more robust with CMs potentially close to their 
	% setpoint limits.
	%
	% INPUTS
	% ------
	% `SC`:: SC base structure
	% `mOrd`:: 	Ordinate of BBA magnet
	% `BPMord`:: Ordinate of BBA BPM 
	% `nDim`:: Hor/Ver plane
	% `par`:: Auxiliary structure
	%
	% RETURN VALUES
	% -------------
	% `CMords`:: CM ordinates for creating orbit bump
	% `CMvec`:: CM setpoint matrix for orbit bumps in BBA measurement

		
	% Exclude dip. compensation CM (needed if considered quadrupole has reverse bending)
	tmpCMind = find(par.RMstruct.CMords{1}==mOrd);
	if ~isempty(tmpCMind)
		par.RMstruct.RM(:,tmpCMind)      = [];
		par.RMstruct.CMords{1}(tmpCMind) = [];
	end
	
	% Get actual index-ordinate pairing of BBA-BPM and orbit feedback BPMs
	tmpBPMind = find(BPMord==par.RMstruct.BPMords);
		
	% Define reference orbit 
	if par.useBPMreadingsForOrbBumpRef
		% Use actual BPM readings and add BBA-BPM target offset
		R0 = SCgetBPMreading(SC);
		R0(nDim,tmpBPMind) = R0(nDim,tmpBPMind) + par.BBABPMtarget;
	else
		% Zero everywhere except at BBA-BPM
		R0 = zeros(2,length(par.RMstruct.BPMords));
		R0(nDim,tmpBPMind) = par.BBABPMtarget;
	end
	
	% Define BPM weighting factors (empirically found to work sufficiently well)
	W0 = ones(2,length(par.RMstruct.BPMords));
	W0(nDim,max(1,tmpBPMind-par.orbBumpWindow):(tmpBPMind-1)) = 0;
	W0(nDim,(tmpBPMind+1):min(length(par.RMstruct.BPMords),tmpBPMind+par.orbBumpWindow)) = 0;
	
	% Run feedback
	[CUR,~] = SCfeedbackRun(SC,par.RMstruct.MinvCO,...
									'weight',[W0(1,:)'; W0(2,:)'],...
									'R0',[R0(1,:)'; R0(2,:)'],...
									'target',0,...
									'maxsteps',50,...
									'scaleDisp',par.RMstruct.scaleDisp,...
									'verbose',par.verbose,...
									'BPMords',par.RMstruct.BPMords,...
									'CMords',par.RMstruct.CMords,...
									'eps',1E-6); 
		
	
	% Read initial and current CM setpoints
	for nDim=1:2
		for nCM=1:length(par.RMstruct.CMords{nDim})
			CMvec0{nDim}(nCM) = SCgetCMSetPoints(SC,par.RMstruct.CMords{nDim}(nCM),nDim);
			deltaCM{nDim}(nCM) = SCgetCMSetPoints(SC,par.RMstruct.CMords{nDim}(nCM),nDim) - SCgetCMSetPoints(CUR,par.RMstruct.CMords{nDim}(nCM),nDim);
		end
	end
		
	% Create matrix of CM setpoints used in BBA measurement
	factor = linspace(-1,1,par.nSteps);
	for nDim=1:2
		for nStep=1:par.nSteps
			CMvec{nDim}(nStep,:) = CMvec0{nDim} + factor(nStep) * deltaCM{nDim};
		end
	end
	% Write CM ordinates for function output
	CMords = par.RMstruct.CMords;
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot BBA step
function plotBBAstep(SC,BPMind,jBPM,nDim,nQ,mOrd,nKick,par)
	% This function plots the particle trajectories/orbit at each measurement step.
	% Note: it is recommended to adjust properties such as plot limits to the users discretion.
	
	global plotFunctionFlag
	% Get s-positions
	sPos = findspos(SC.RING,1:length(SC.RING))';
	% Axes limits
	% 	xLim = [0 12.5];yLim = .3*[-1 1];
	xLim = sPos(mOrd)+[-10 10];yLim = 1.3*[-1 1];
	
	% Clear figure at beginning of measurement
	if nQ==1 && nKick==1;figure(99);clf;end
	% Get trajectories
	plotFunctionFlag=1;[B,T]=SCgetBPMreading(SC);plotFunctionFlag=[];
	% Select only first particle
	T=SCparticlesIn3D(T,SC.INJ.nParticles);T=T(:,:,1);
	% Select axes
	figure(99);subplot(length(par.magSPvec{nDim,jBPM}),1,nQ);hold on;
	% Plot downstream BPM readings
	plot(sPos(SC.ORD.BPM),1E3*B(nDim,1:length(SC.ORD.BPM)),'o');
	% PLot BBA-BPM readings
	plot(sPos(SC.ORD.BPM(BPMind)),1E3*B(nDim,BPMind),'ko','MarkerSize',10,'MarkerFaceColor','k');
	% Plot trajectories
	plot(sPos,1E3*T(2*nDim-1,1:length(SC.RING)),'-');
	% Plot BBA magnet
	rectangle('Position',[sPos(mOrd),-1,sPos(mOrd+1)-sPos(mOrd),1 ],'FaceColor',[0 0.4470 0.7410]);
	xlim(xLim);ylim(yLim)
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot BBA results
function plotBBAResults(SC,initOffsetErrors,errorFlags,jBPM,BPMords,magOrds)
	% Plots the initial and final BPM offsets with respect to the magnet centers
	
	
	
	% Define figure of merit
	fom0  = initOffsetErrors;
	fom = 	getBPMoffsetFromMag(SC,BPMords,magOrds);
	
	% Remove BPMs which have not yet been measured
	fom(:,jBPM+1:end) = nan;
	
	if size(BPMords,2)==1
		nSteps = 1;
	else
		nSteps = 1.1*max(abs(fom0(:)))*linspace(-1,1,floor(size(BPMords,2)/3));
	end
	
	figure(90),clf;tmpCol=get(gca,'colororder');
	subplot(3,1,1);hold on
	
	% Plot histogram of new BPM offsets
	for nDim=1:size(BPMords,1)
		[a,b] = hist(fom(nDim,:),nSteps);
		plot(1E6*b,a,'LineWidth',2);
	end
	% Plot histogram of initial BPM offsets
	[a,b] = hist(fom0(:),nSteps);
	plot(1E6*b,a,'k-','LineWidth',2)
	
	if size(BPMords,1)>1
		legend({sprintf('Horizontal rms: $%.0f\\mu m$',1E6*sqrt(mean(fom(1,:).^2,'omitnan'))),sprintf('Vertical rms: $%.0f\\mu m$',1E6*sqrt(mean(fom(2,:).^2,'omitnan'))),sprintf('Initial rms: $%.0f\\mu m$',1E6*sqrt(mean(fom0(:).^2,'omitnan')))},'Interpreter','latex')
	end
	xlabel('Final BPM offset w.r.t. magnet [$\mu$m]'),ylabel('Number of counts');set(gca,'box','on')

	
	subplot(3,1,2);hold on;p=zeros(1,4);

	% Plot new BPM offset w.r.t. magnet centers
	for nDim=1:size(BPMords,1)
		x = find(ismember(SC.ORD.BPM,BPMords(nDim,:)));
		if any(~isnan(fom(nDim,errorFlags(nDim,:)==0)))
			p(nDim)=plot(x(errorFlags(nDim,:)==0),1E6*abs(fom(nDim,errorFlags(nDim,:)==0)),'O','LineWidth',2,'Color',tmpCol(nDim,:));
		end
		if any(~isnan(fom(nDim,errorFlags(nDim,:)~=0)))
			p(2+nDim)=plot(x(errorFlags(nDim,:)~=0),1E6*abs(fom(nDim,errorFlags(nDim,:)~=0)),'X','LineWidth',2,'Color',tmpCol(nDim,:));
		end
	end
	ylabel('Final offset [$\mu$m]'),xlabel('Index of BPM');set(gca,'XLim',[1 length(SC.ORD.BPM)],'box','on')
	legStr = {'Horizontal','Vertical','Horizontal failed','Vertical failed'};
	legend(p(p~=0),legStr{p~=0});

	
	subplot(3,1,3);hold on;p=zeros(1,4);

	% Plot offset change
	x = find(ismember(SC.ORD.BPM,BPMords(nDim,:)));
	for nDim=1:size(BPMords,1)
		plot(x,1E6*(fom0(nDim,:)-fom(nDim,:)),'d','LineWidth',2)
	end
	ylabel('Offsets change [$\mu$m]'),xlabel('Index of BPM');set(gca,'XLim',[1 length(SC.ORD.BPM)],'box','on')
	legend({'Horizontal','Vertical'});

	set(gcf,'color','w');set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');set(findall(gcf,'-property','FontSize'),'FontSize',18);
	drawnow
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get BPM offset with respect to specified magnet centers
function offset = getBPMoffsetFromMag(SC,BPMords,magOrds)
	% Get BPM offsets with respect to magnet centers
	%
	% INPUTS
	% ------
	% `SC`:: SC base structure
	% `BPMords`:: BPM ordinates used in BBA
	% `magOrds`:: Magnet ordinates used in BBA
	%
	% RETURN VALUES
	% -------------
	% `offset`:: BPM offsets with respect to magnet centers
	
	% Allocate output
	offset = nan(size(BPMords));
	
	% Loop over BPMs in both planes
	for nDim=1:size(BPMords,1)
		for nBPM=1:size(BPMords,2)
			offset(nDim,nBPM) = SC.RING{BPMords(nDim,nBPM)}.Offset(nDim) + SC.RING{BPMords(nDim,nBPM)}.SupportOffset(nDim) - SC.RING{magOrds(nDim,nBPM)}.MagnetOffset(nDim) - SC.RING{magOrds(nDim,nBPM)}.SupportOffset(nDim);
		end
	end

end





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fake measurement to remove outliers
function SC = fakeMeasurement(SC,BPMords,magOrds,errorFlags)
	% This function intends to mimic the operater's ability to identify errors in the measurement procedure,
	% adjust the fine tuning parameters for individual BPMs and repeat the measurement.
	% Here, the rms value of the difference between the BPM offsets and the magnet centers is calculated 
	% for both planes for the succesfull BPMs.
	% All BPM offsets at which the BBA procedure failed are removed and artificially
	% re-generated using a Gaussian distribution with two sigma cutoff with an rms value as described abive.
	%
	% INPUTS
	% ------
	% `SC`:: SC base structure
	% `BPMords`:: BPM ordinates used in BBA
	% `magOrds`:: Magnet ordinates used in BBA
	% `errorFlags`:: Error flags for the BBA measruement
	%
	% RETURN VALUES
	% -------------
	% `SC`:: SC base structure with removed outliers

	
	% Get BPM offset errors w.r.t. magnet centers
	finOffsetErrors = getBPMoffsetFromMag(SC,BPMords,magOrds);

	% Remove failed BPMs
	finOffsetErrors(errorFlags~=0) = nan;

	% PRINT INFORMATION
	fprintf('Final offset error is %.1f|%.1f um (hor|ver) with %d|%d measurement failures -> being re-calculated now.\n' ,1E6*sqrt(mean(finOffsetErrors.^2,2,'omitnan')),sum(errorFlags~=0,2))

	% Loop over BPMs in both planes
	for nBPM=1:size(BPMords,2)
		for nDim=1:2
			% Check if measurement failed
			if errorFlags(nDim,nBPM)~=0
				% Calculate fake measurement result
				fakeBPMoffset = SC.RING{magOrds(nDim,nBPM)}.MagnetOffset(nDim) + SC.RING{magOrds(nDim,nBPM)}.SupportOffset(nDim) - SC.RING{BPMords(nDim,nBPM)}.SupportOffset(nDim) + sqrt(mean(finOffsetErrors(nDim,:).^2,'omitnan')) * SCrandnc(2);
				% Ensure no NaN
				if ~isnan(fakeBPMoffset)
					% Write new BPM offset
					SC.RING{BPMords(nDim,nBPM)}.Offset(nDim) = fakeBPMoffset;
				else
					fprint('BPM offset not reasigned, NaN.\n')
				end
			end
		end
	end
end

