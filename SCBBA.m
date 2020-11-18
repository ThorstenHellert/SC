function [SC,BPMstd] = SCBBA(SC,BPMords,QuadOrds,varargin)
% SCBBA
% =====
%
% NAME
% ----
% SCBBA - Performs model independend beam based alignment
%
% SYNOPSIS
% --------
% `[SC, BPMstd] = SCBBA(SC, BPMords, QuadOrds, [, options])`
%
%
% DESCRIPTION
% -----------
% Perform a model independend beam based alignment procedure using either two-turn
% trajectories or closed orbit bumps.
% In two-turn mode, for each BPM the injected beam trajectory is varied within a
% user defined range and the corresponding quadrupole is exercised on a user
% defined vector of setpoints. If the initial trajectory variation causes beam
% losses, the trajectory variation is reduced until all injected beam
% trajectories reach the BBA-BPM. If the final trajectory variation at the
% considered BBA-BPM is below a user defined threshold, a quadrupole may be
% exercised to change the phase advance between the injection point and the
% BPM. 
% In orbit mode an orbit bump is generated at each considerded BPM using orbit feedback with weighting factors
% Finally, for each injected beam trajectory or orbit bump step the offset variation due
% the change in quadrupole strength is recorded at the BPMs used for measurement. A line
% fit is used to determine the center of the BBA-BPM with respect to the used
% quadrupole.
%
% INPUTS
% ------
%
% `SC`::
%	SC base structure
% `BPMords`::
%	[2 x n] array of BPM ordinates
% `QuadOrds`::
%	[2 x n] array of quadrupole ordinates for the corresponding BPMs in `BPMords`
%
% RETURN VALUES
% -------------
%
% `SC`::
%	SC base structure with updated BPM offsets
% `BPMstd`::
%	[2 x n] array of standard deviations of the fitted BPM centers
%		
%  
% OPTIONS
% -------
% The following options can be specified as name-value pairs:
%
% `'mode'` (`'ORB'`)::
%   Orbit or two-turn trajectory mode (`'2turn'`)
% `'removeOutliersAt'` (`Inf`)::
%   If the calculated BPM offset w.r.t. the magnet center is above the specified value, the measruement is discarted and the BPM offset is not updated.
% `'removeOutliersAtSigma'` (`Inf`)::
%   This option intends to mimic the operater's ability to identify errors in the measurment procedure and adjust the fine tuning parameters for individual BPMs.
%	After performing the measurement routine, the rms value of the difference between the BPM offsets and the magnet centers is calculated for both planes.
%   All offset errors outside the number of sigmas specified here are removed and artificially generated using a Gaussian distribution with two sigma cutoff with the rms value as described abive.
% `'BPMrangeAtBBABBPMTarget'` (`1E-3`)::
%	Required offset variation at BBA-BPM (BPM adjacent to quad)
%	which should be achieved by trajectory change or orbit bump.
% `'minBPMrangeAtBBABBPM'` (`0.5E-3`)::
%	Threshold of change of offset at BBA-BPM, if below BBA
%	evaluation is not performed.
% `'minBPMrangeDownstream'` (`100E-6`)::
%	Minimum change of offset at downstream BPMs to be included
%	in calculation.
% `'maxKickatInjection'` (`0.9E-3`)::
%	 Maximum kick [rad] change at injection (2-turn mode only).
% `'maxOffsetatInjection'` (`0.9E-3`)::
%	 Maximum offset [m] change at injection (2-turn mode only).
% `'nSteps'` (`10`)::
%	 Number of different trajectories/orbit bump steps for each quadrupole setting
% `'quadSPvec'` (`[0.95 1.05]`)::
%	 Strength variation of quadrupoles.
% `'quadSPflag'` (`'rel'`)::
%	 Specify how quadrupole setpoint should be changed, e.g. 'rel' or 'abs'.
% `'switchOffSext'` (`0`)::
%	 Flag specifying if sextupole coil in BBA magnet should be switched off 
%    (e.g. if quadrupole trim coils are used).
% `'RMstruct'` (`[]`)::
%	 Structure containing pre-calculated response matrix etc. for orbit feedback (orbit mode only).
%    If empty the relevant arrays are calculated with default parameters.
% `'quadOrdPhaseAdvance'` (`[]`)::
%	 Quadrupole ordinate which is used to change the phase advance between
%	 injection and BBA-BPM (2-turn mode only).
% `'quadStrengthPhaseAdvance'` (`[0.95 0.8 1.05]`)::
%	 Relative quadrupole strength variation used to change the phase
%	 advance between injection and BBA-BPM (2-turn mode only).
% `'nXPointsNeededAtMeasBPM'` (`3`)::
%	 Number of x-positions at BBA-BPM required for linear regression at
%	 measurement BPMs.
% `'maxNumOfDownstreamBPMs'` (`length(SC.ORD.BPM)`)::
%	 Number downstream BPMs considered in the data evaluation (2-turn mode only).
% `'minSlopeForFit'` (`0.03`)::
%	 Minimum fitted slope at measurement BPMs with respect to quadrupole
%	 change which is still used in determining BPM center (a very small slope 
%    usually indicates an unfortunate phase advance and can significantly affect 
%    the measurement if not excluded). 
% `'maxStdForFittedCenters'` (`600E-6`)::
%	 If standard deviation of the fitted BPM centers as determined by all
%	 downstream BPMs exceeds this value, output will be `'nan'`.
% `'plotLines'` (0)::
%	If true, each injected beam and intermediate BBA results will be
%	plotted.
% `'plotResults'` (0)::
%	If true, final BBA results are plotted.
% `'verbose'` (0)::
%	If true, debug information is printed.
%
%
% SEE ALSO
% --------
% *SCsetMags2SetPoints*, *SCgetBPMreading*, *SCfeedbackRun*

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Input check
	
	p = inputParser;
	addOptional(p,'mode','ORB');
	addOptional(p,'removeOutliersAt',Inf);
	addOptional(p,'removeOutliersAtSigma',Inf);
	addOptional(p,'BPMrangeAtBBABBPMTarget',1E-3);
	addOptional(p,'minBPMrangeAtBBABBPM',500E-6);
	addOptional(p,'minBPMrangeOtherBPM',100E-6);
	addOptional(p,'maxTrajChangeAtInjection',[.9E-3 .9E-3]);
	addOptional(p,'nSteps',10);
	addOptional(p,'quadSPvec',[0.95,1.05]);
	addOptional(p,'quadSPflag','rel');
	addOptional(p,'switchOffSext',0);
	addOptional(p,'RMstruct',[]);
	addOptional(p,'quadOrdPhaseAdvance',[ ]);
	addOptional(p,'quadStrengthPhaseAdvance',[0.95 0.8 1.05]);
	addOptional(p,'nXPointsNeededAtMeasBPM',3);
	addOptional(p,'maxNumOfDownstreamBPMs',length(SC.ORD.BPM));
	addOptional(p,'minSlopeForFit',0.03);
	addOptional(p,'maxStdForFittedCenters',600E-6);
	addOptional(p,'plotLines',0);
	addOptional(p,'plotResults',0);
	addOptional(p,'verbose',0);
	parse(p,varargin{:});
	par = p.Results;

	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Initialization

	% Save initial BPM offset errors (for plotting)
	initOffsetErrors = getBPMoffsetFromMag(SC,BPMords,QuadOrds);
	
	% Offset and kick angle variation at injection (for trajectory mode)
	kickVec0  = par.maxTrajChangeAtInjection' .* repmat(linspace(-1,1,par.nSteps),2,1);

	% Save initial injection trajectory (for trajectory mode)
	initialZ0 = SC.INJ.Z0;
	
	% Calculateall relevant arrays for orbit correction if not supplied by user
	if strcmp(par.mode,'ORB')
		if isempty(par.RMstruct)
			fprintf('Calcualting orbit respone matrix and dispersion.\n')
			par.RMstruct.RM        = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'trackMode','ORB','useIdealRing',1); 
			par.RMstruct.eta       = SCgetModelDispersion(SC,SC.ORD.BPM,SC.ORD.Cavity);
			par.RMstruct.scaleDisp = 1E7;
			par.RMstruct.BPMords   = SC.ORD.BPM;
			par.RMstruct.CMords    = SC.ORD.CM;
			par.RMstruct.alpha     = 1;
		end
	end
	
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Perform BBA routine
	
	% Loop over BPMs
	for jBPM=1:size(BPMords,2) % jBPM: Index of BPM adjacent to quad for BBA

		% Horizontal/vertical
		for nDim=1:size(BPMords,1)
			if par.verbose;fprintf('BBA-BPM %d/%d, nDim = %d',jBPM,size(BPMords,2),nDim);end

			% Save initial machine state (for convenience)
			SC0 = SC;
			
			% Get BPM index-ordinate pairing
			BPMind = find(BPMords(nDim,jBPM)==SC.ORD.BPM);
			
			% Define ordinate of quadrupole for BBA
			qOrd = QuadOrds(nDim,jBPM);


			% Switch off sextupole coil at BBA magnet?
			if par.switchOffSext
				% Switch off sextupole coil
				SC = SCsetMags2SetPoints(SC,qOrd,2,3,0,'method','abs');
				% Get pseudo inverse
				MinvCO = SCgetPinv([par.RMstruct.RM par.RMstruct.scaleDisp*par.RMstruct.eta],'alpha',par.RMstruct.alpha);
				% Run orbit feedback
				[SC,~] = dev_SCfeedbackRun(SC,MinvCO,...
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
					[CMords,CMvec] = getOrbitBump(SC,qOrd,BPMords(nDim,jBPM),nDim,par);
					
					% Perform data measrurment
					[BPMpos,tmpTra] = dataMeasurement(SC,qOrd,BPMind,nDim,par,CMords,CMvec);

				case '2turn'
					% Scale maximum initial trajectory variation before beam gets lost
					[kickVec, BPMrange] = scaleInjectionToReachBPM(SC,BPMind,nDim,initialZ0,kickVec0,par);
					
					% Check if offset range at BBA-BPM is to small and vary phase advance
					if ~isempty(par.quadOrdPhaseAdvance) && BPMrange < par.BPMrangeAtBBABBPMTarget
						[SC,kickVec] = scanPhaseAdvance(SC,BPMind,nDim,initialZ0,kickVec0,par);
					end
					
					% Perform data measurement
					[BPMpos,tmpTra] = dataMeasurement(SC,qOrd,BPMind,nDim,par,initialZ0,kickVec);
			end
			
			% Perform data evaluation
			[BPMcenter,BPMstd] = dataEvaluation(SC,BPMords,jBPM,BPMpos,tmpTra,nDim,qOrd,par);

			% Reset machine (should ideally be done by setting all magnet setpoints to intial)
			SC = SC0;
	
			% Outlier detection
			if  SC.RING{BPMords(nDim,jBPM)}.Offset(nDim) + BPMcenter - SC.RING{qOrd}.MagnetOffset(nDim) > par.removeOutliersAt
				BPMcenter = nan;
			end
			
			% Check if BPM center could be identified in measurement
			if ~isnan(BPMcenter)
				% Apply BBA results
				SC.RING{BPMords(nDim,jBPM)}.Offset(nDim) = SC.RING{BPMords(nDim,jBPM)}.Offset(nDim) + BPMcenter;
			end
		end

		% Plot BBA results
		if par.plotResults
			plotBBAResults(SC,initOffsetErrors,jBPM,BPMords,QuadOrds)
		end
	end
	
	% Remove outliers
	SC = removeOutliers(SC,BPMords,QuadOrds,par);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data measurement
function [BPMpos,tmpTra] = dataMeasurement(SC,qOrd,BPMind,nDim,par,varargin)
	% This function performs the data measruement by looping through the quadrupole
	% setpoints and changing either the injected beam trajectory or the CMs to create 
	% a local orbit bump.
	%
	% INPUTS
	% ------
	% `SC`:: SC base structure
	% `qOrd`:: 	ordinate of BBA quadrupole
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

	switch par.mode
		case 'ORB'
			
			% Used CM ordinates and setpoint matrix
			CMords = varargin{1};
			CMvec  = varargin{2};
			
			% Measurment steps per quad setting
			nMsteps = size(CMvec{nDim},1);

			% Prealocate x and y data
			tmpTra = nan(nMsteps,length(par.quadSPvec),length(SC.ORD.BPM));
			BPMpos = nan(nMsteps,length(par.quadSPvec));

		case '2turn'
			% Initial trajectory and kick variations
			initialZ0 = varargin{1};
			kickVec   = varargin{2};
			
			% Measurment steps per quad setting
			nMsteps = size(kickVec,2);
			
			% Prealocate x and y data
			tmpTra = nan(nMsteps,length(par.quadSPvec),par.maxNumOfDownstreamBPMs);
			BPMpos = nan(nMsteps,length(par.quadSPvec));
	end
	
	% Loop over different quadrupole settings
	for nQ=1:length(par.quadSPvec)

		% Set quadrupole to different setpoints
		SC = SCsetMags2SetPoints(SC,qOrd,2,2,par.quadSPvec(nQ),...
			'method',par.quadSPflag,...
			'dipCompensation',1);

		% Loop over different trajectories
		for nKick=1:nMsteps

			switch par.mode
				case 'ORB'
					% Set CM
					for nD=1:2
						SC = SCsetCMs2SetPoints(SC,CMords{nD},CMvec{nD}(nKick,:),nD,'abs');
					end
				case '2turn'
					% Vary initial trajectory
					SC.INJ.Z0(2*nDim)   = initialZ0(2*nDim  ) + kickVec(2,nKick); % kick angle
					SC.INJ.Z0(2*nDim-1) = initialZ0(2*nDim-1) + kickVec(1,nKick); % offset
			end
			
			% Calculate beam reading
			B = SCgetBPMreading(SC);

			% JUST FOR PLOTTING
			if par.plotLines; plotBBAstep(SC,BPMind,nDim,nQ,qOrd,nKick,par);end

			% Save BBA-BPM readings
			BPMpos(nKick,nQ) = B(nDim,BPMind);
	
			switch par.mode
				case 'ORB'
					% Save all BPM readings
					tmpTra(nKick,nQ,:) = B(nDim, : );
				case '2turn'
					% Save downstream BPM readings
					tmpTra(nKick,nQ,:) = B(nDim, (BPMind+1):(BPMind+par.maxNumOfDownstreamBPMs) );
			end
		end
	end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data evaluation
function [BPMcenter,BPMstd] = dataEvaluation(SC,BPMords,jBPM,BPMpos,tmpTra,nDim,qOrd,par)
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
	% `BPMcenter`:: Calculated BPM offset (to be on lattice element)
	% `BPMstd`:: Calculated BPM offset error 

	
	
	% Plot measurement details
	if par.plotLines;figure(56),clf;hold on; p(1)=plot3(0,1E6*SC.RING{qOrd}.T2(2*nDim-1),0,'rD','MarkerSize',40,'MarkerFaceColor','b'); hold on;set(gca,'box','on');end
 
	tmpCenter = nan(1,(size(tmpTra,2)-1)*par.maxNumOfDownstreamBPMs);
	tmpRangeX = zeros(1,(size(tmpTra,2)-1)*par.maxNumOfDownstreamBPMs);
	tmpRangeY = zeros(1,(size(tmpTra,2)-1)*par.maxNumOfDownstreamBPMs);

	i = 0;
	% Loop over downstream BPMs
	for nBPM=1:par.maxNumOfDownstreamBPMs

		% y-data: downstream trajectory differences with respect to quad setting
		y0 = diff(tmpTra(:,:,nBPM),1,2);
		% x data: BBA-BPM readings
		x0 = repmat(mean(BPMpos,2),1,size(y0,2));

		% Loop over different quadrupole differences
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

			% Check if sufficient transmission and signal is achieved
			if length(x)>=par.nXPointsNeededAtMeasBPM ...        % More than 3 different points at BBA-BPM
					&& tmpRangeX(i) > par.minBPMrangeAtBBABBPM ...     % Range at BBA-BPM        > par.minBPMrangeAtBBABBPM
					&& tmpRangeY(i) > par.minBPMrangeOtherBPM        % Range at downstream BPM > pra.minBPMrangeOtherBPM

				% Linear regression
				sol = [ones(size(x)),x]\y;
				offset = sol(1);
				slope  = sol(2);

				% Check if fitted slope is to small
				if abs(slope) < par.minSlopeForFit
					slope = nan;
				end

				% Calculate fitted magnetic center
				tmpCenter(i) = -offset/slope;
			else
				tmpCenter(i) = nan;
				slope=nan;offset=nan; % just needed for plotting
			end

			% Just plotting
			if par.plotLines 
				p(2)=plot3(repmat(jBPM+nBPM,size(x)),1E6*x,1E3*y,'ko');
				tmp=plot3(repmat(jBPM+nBPM,size(BPMpos)),1E6*BPMpos,1E3*(slope*BPMpos+offset),'k-');p(3)=tmp(1);
				p(4)=plot3(jBPM+nBPM,1E6*tmpCenter(nBPM),0,'Or','MarkerSize',10);
			end

		end
	end

	% Check if measurment requirements have been met
	if (max(tmpRangeX) < par.minBPMrangeAtBBABBPM) ...                % Max. range at BBA-BPM to small
			|| (max(tmpRangeY) < par.minBPMrangeOtherBPM) ...       % Max. range at downstream BPM to small
			|| std(tmpCenter,'omitnan') > par.maxStdForFittedCenters  % Fitted magnetic centers to far spread out

		BPMcenter = nan;
		BPMstd    = nan;
	else
		% Calculate mean magnetic center and error
		BPMcenter = mean(tmpCenter,'omitnan');
		BPMstd    = std( tmpCenter,'omitnan');
	end

	% Just plotting
	if par.plotLines
		p(5)=plot3(0,1E6*BPMcenter,0,'kD','MarkerSize',30,'MarkerFaceColor','r'); 
		p(6)=plot3(0,1E6*(SC.RING{BPMords(nDim,jBPM)}.Offset(nDim)+BPMcenter),0,'kD','MarkerSize',30,'MarkerFaceColor','g'); 
		title(sprintf('BBA-BPM: %d // qOrd: %d // qFam: %s // nDim: %d // FinOffset = %.0f $\\mu m$',jBPM,qOrd,SC.RING{qOrd}.FamName,nDim,1E6*abs(SC.RING{BPMords(nDim,jBPM)}.Offset(nDim) + SC.RING{BPMords(nDim,jBPM)}.SupportOffset(nDim) + BPMcenter - SC.RING{qOrd}.MagnetOffset(nDim) - SC.RING{qOrd}.SupportOffset(nDim))));
		legend(p,{'Magnet center','Measured offset change','Line fit','Fitted BPM offset (individual)','Fitted BPM offset (final)','Corrected BPM offset'})
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
	qOrd = par.quadOrdPhaseAdvance;
	qVec = par.quadStrengthPhaseAdvance;
	q0   = SC.RING{qOrd}.SetPointB(2);

	% Initialize all BPM ranges
	allBPMRange = zeros(length(qVec));

	% Loop over quadrupole setpoints
	for nQ=1:length(qVec)
		if par.verbose;fprintf('BBA-BPM range to small, try to change phase advance with quad ord %d to %.2f of nom. SP.\n',par.quadOrdPhaseAdvance,qVec(nQ));end

		% Set quadrupole to different setpoints
		SC = SCsetMags2SetPoints(SC,qOrd,2,2,qVec(nQ),'method','rel','dipCompensation',1);

		% Rescale initial trajectory variation
		[kickVec, BPMrange] = scaleInjectionToReachBPM(SC,BPMind,nDim,initialZ0,kickVec0,par);

		% Record BPM range
		allBPMRange(nQ) = BPMrange;

		if par.verbose;fprintf('Initial trajectory variation scaled to [%.0f|%.0f]%% of its initial value, BBA-BPM range %.0fum.\n',100*(kickVec([1 end])./kickVec0([1 end])),1E6*BPMrange);end

		% End loop if offset range at BBA-BPM is sufficient
		if ~( BPMrange < par.BPMrangeAtBBABBPMTarget )
			break
		end
	end

	% Check if offset range at BBA-BPM is still to small
	if BPMrange < par.BPMrangeAtBBABBPMTarget
		% Check if any BPM range is better than the initial one
		if BPMrange<max(allBPMRange)
			[~,nBest] = max(allBPMRange);
			% Set quadrupole to best setpoints
			SC = SCsetMags2SetPoints(SC,qOrd,2,2,qVec(nBest),'method','rel','dipCompensation',1);
			if par.verbose;fprintf('Changing phase advance of quad with ord %d NOT succesfull, returning to best value with BBA-BPM range = %.0fum.\n',qOrd,1E6*max(allBPMRange));end
		else
			% Reset quadrupole to initial setpoints
			SC = SCsetMags2SetPoints(SC,qOrd,2,2,q0,'method','abs','dipCompensation',1);
			if par.verbose;fprintf('Changing phase advance of quad with ord %d NOT succesfull, returning to initial setpoint.\n',qOrd);end
		end
	else
		if par.verbose;fprintf('Change phase advance of quad with ord %d successful. BBA-BPM range = %.0fum.\n',qOrd,1E6*BPMrange);end
	end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get orbit bump
function [CMords,CMvec] = getOrbitBump(SC,qOrd,BPMord,nDim,par)
	% This function is used to create a pseudo orbit bump by running orbit feedback on the actual 
	% machine with target 
	% of zero except at the BBA-BPM (here the target is 'par.BPMrangeAtBBABBPMTarget') and a 
	% few BPMs upstream and downstream (note that this can be 
	% done better by creating a 'real' orbut bump using 3 CMs, however we found that this 
	% procedure is faster and more robust with limited CM strenghts.
	%
	% INPUTS
	% ------
	% `SC`:: SC base structure
	% `qOrd`:: 	ordinate of BBA quadrupole
	% `BPMord`:: Ordinate of BBA BPM 
	% `nDim`:: Hor/Ver plane
	% `par`:: Auxiliary structure
	%
	% RETURN VALUES
	% -------------
	% `CMords`:: CM ordinates for creating orbit bump
	% `CMvec`:: CM setpoint matrix for orbit bumps in BBA measurement

		
	% Exclude dip. compensation CM (needed if consideed quadrupole has reverse bending)
	par.RMstruct.CMords{1}(find(par.RMstruct.CMords{1}==qOrd)) = [];
	
	% Get actual index-ordinate pairing of BBA BPM and orbit feedback BPMs
	tmpBPMind = find(BPMord==par.RMstruct.BPMords);
		
	% Define reference orbit (zero everywhere except at BBA BPM)
	R0 = zeros(2,length(par.RMstruct.BPMords));
	R0(nDim,tmpBPMind) = par.BPMrangeAtBBABBPMTarget;
	
	% Define BPM weighting factors (empirically found to work sufficiently well)
	window = 5;
	W0 = ones(2,length(par.RMstruct.BPMords));
	W0(nDim,max(1,tmpBPMind-window):(tmpBPMind-1)) = 0;
	W0(nDim,(tmpBPMind+1):min(length(par.RMstruct.BPMords),tmpBPMind+window)) = 0;
	
	
	% Get pseudo inverse
	MinvCO = SCgetPinv([par.RMstruct.RM par.RMstruct.scaleDisp*par.RMstruct.eta],'alpha',par.RMstruct.alpha);
	% Run feedback
	[CUR,~] = SCfeedbackRun(SC,MinvCO,...
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
	factor = 1*linspace(-1,1,par.nSteps);
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
function plotBBAstep(SC,BPMind,nDim,nQ,qOrd,nKick,par)
	% This function plots the particle trajectories/orbit at each measurement step.
	% Note: it is recommended to adjust properties such as plot limits to the users discretion.
	
	global plotFunctionFlag
	xLim = [0 12.5];yLim = 3*[-1 1];
	% Clear figure at beginning of measurment
	if nQ==1 && nKick==1;figure(99);clf;end
	% Get s-positions
	sPos = findspos(SC.RING,1:length(SC.RING))';
	% Get trajectories
	plotFunctionFlag=1;[B,T]=SCgetBPMreading(SC);plotFunctionFlag=[];
	% Select only first particle
	T=SCparticlesIn3D(T,SC.INJ.nParticles);T=T(:,:,1);
	% Select axes
	figure(99);subplot(length(par.quadSPvec),1,nQ);hold on;
	% Plot downstream BPM readings
	plot(sPos(SC.ORD.BPM),1E3*B(nDim,1:length(SC.ORD.BPM)),'o');
	% PLot BBA-BPM readings
	plot(sPos(SC.ORD.BPM(BPMind)),1E3*B(nDim,BPMind),'ko','MarkerSize',10,'MarkerFaceColor','k');
	% Plot trajectories
	plot(sPos,1E3*T(2*nDim-1,1:length(SC.RING)),'-');
	% Plot quadrupole magnet
	rectangle('Position',[sPos(qOrd),-1,sPos(qOrd+1)-sPos(qOrd),1 ],'FaceColor',[0 0.4470 0.7410]);
	xlim(xLim);ylim(yLim)
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot BBA results
function plotBBAResults(SC,initOffsetErrors,jBPM,BPMords,QuadOrds)
	% Plots the initial and final BPM offsets with respect to the magnet centers
	
	
	
	% Define figure of merit
	fom0  = initOffsetErrors;
	fom = 	getBPMoffsetFromMag(SC,BPMords,QuadOrds);
	
	% Remove BPMs which have not yet been measured
	fom(:,jBPM+1:end) = nan;
	
	if size(BPMords,2)==1
		nSteps = 1;
	else
		nSteps = 1.1*max(abs(fom0(:)))*linspace(-1,1,floor(size(BPMords,2)/3));
	end
	
	figure(90),clf;
	subplot(2,1,1);hold on
	
	% Plot histogram of new BPM offsets
	for nDim=1:size(BPMords,1)
		[a,b] = hist(fom(nDim,:),nSteps);
		plot(1E6*b,a,'LineWidth',2);
	end
	% Plot histogram of initial BPM offsets
	[a,b] = hist(fom0(:),nSteps);
	plot(1E6*b,a,'k-','LineWidth',2)
	
	if size(BPMords,1)>1
		legend({sprintf('Horizontal rms: $%.0f\\mu m$',1E6*sqrt(mean(fom(1,:).^2,'omitnan'))),sprintf('Vertical rms: $%.0f\\mu m$',1E6*sqrt(mean(fom(2,:).^2,'omitnan'))),sprintf('Uncorrected rms: $%.0f\\mu m$',1E6*sqrt(mean(fom0(:).^2,'omitnan')))},'Interpreter','latex')
	end
	xlabel('BPM offsets w.r.t. quads [$\mu$m]'),ylabel('Number of counts');set(gca,'box','on')

	subplot(2,1,2);hold on

	% Plot new BPM offset w.r.t. magnet centers
	for nDim=1:size(BPMords,1)
		x = find(ismember(SC.ORD.BPM,BPMords(nDim,:)));
		plot(x,1E6*abs(fom(nDim,:)),'O','LineWidth',2)
	end
	ylabel('BPM offsets w.r.t. quads [$\mu$m]'),xlabel('Index of BPM');set(gca,'XLim',[1 length(SC.ORD.BPM)],'box','on')


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
	% `qOrd`:: 	ordinate of BBA quadrupole
	% `BPMords`:: BPM ordinates used in BBA
	% `QuadOrds`:: Quadrupole ordinates used in BBA
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
% Remove outliers
function SC = removeOutliers(SC,BPMords,QuadOrds,par)
	% This function intends to mimic the operater's ability to identify errors in the measurment procedure,
	% adjust the fine tuning parameters for individual BPMs and repeat the measurement.
	% Here, the rms value of the difference between the BPM offsets and the magnet centers is calculated for both planes.
	% All offset errors outside the number of sigmas specified in 'par.removeOutliersAtSigma' are removed and artificially
	% re-generated using a Gaussian distribution with two sigma cutoff with an rms value as described abive.
	%
	% INPUTS
	% ------
	% `SC`:: SC base structure
	% `qOrd`:: 	ordinate of BBA quadrupole
	% `BPMords`:: BPM ordinates used in BBA
	% `QuadOrds`:: Quadrupole ordinates used in BBA
	% `par`:: Auxiliary structure
	%
	% RETURN VALUES
	% -------------
	% `SC`:: SC base structure with removed outliers

	if par.removeOutliersAtSigma==Inf
		return
	end
	
	% Get BPM offset errors w.r.t. magnet centers
	finOffsetErrors = getBPMoffsetFromMag(SC,BPMords,QuadOrds);
	
	% PRINT INFORMATION
	fprintf('Final offset error is %.1f|%.1f um (hor|ver) with %d|%d outliers outside %d-sigma interval -> being removed now.\n' ,1E6*sqrt(mean(finOffsetErrors.^2,2)),sum(abs(finOffsetErrors) > par.removeOutliersAtSigma * sqrt(mean(finOffsetErrors.^2,2)),2),par.removeOutliersAtSigma)
	
	% Loop over BPMs in both planes
	for nBPM=1:size(BPMords,2)
		for nDim=1:2
			% Check if measurement is outside the user defined interval
			if abs(finOffsetErrors(nDim,nBPM)) > par.removeOutliersAtSigma * sqrt(mean(finOffsetErrors(nDim,:).^2))
				% Write new BPM offset 
				SC.RING{BPMords(nDim,nBPM)}.Offset(nDim) = SC.RING{QuadOrds(nDim,nBPM)}.MagnetOffset(nDim) + SC.RING{QuadOrds(nDim,nBPM)}.SupportOffset(nDim) - SC.RING{BPMords(nDim,nBPM)}.SupportOffset(nDim) + sqrt(mean(finOffsetErrors(nDim,:).^2)) * SCrandnc(2);
			end
		end
	end
end	
