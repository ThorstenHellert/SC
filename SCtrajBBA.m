function [SC,BPMstd] = SCtrajBBA(SC,BPMords,QuadOrds,varargin)
% SCBBA
% =====
%
% NAME
% ----
% SCBBA - Performs model independend 2-turn beam based alignment
%
% SYNOPSIS
% --------
% `[SC, BPMstd] = SCBBA(SC, BPMords, QuadOrds, [, options])`
%
%
% DESCRIPTION
% -----------
% Perform a model independend beam based alignment procedure using two-turn
% trajectories. For each BPM the injected beam trajectory is varied within a
% user defined range and the corresponding quadrupole is exercised on a user
% defined vector of setpoints. If the initial trajectory variation causes beam
% losses, the trajectory variation is reduced until all injected beam
% trajectories reach the BBA-BPM. If the final trajectory variation at the
% considered BBA-BPM is below a user defined threshold, a quadrupole may be
% exercised to change the phase advance between the injection point and the
% BPM. Finally, for each injected beam trajectory the trajectory variation due
% the change in quadrupole strength is recorded at the downstream BPMs. A line
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
% `'minBPMrangeAtBBABBPMTarget'` (`1E-3`)::
%	Minimum change of trajectory offset at BBA-BPM (BPM adjacent to quad)
%	which should be achieved by trajectory change.
% `'minBPMrangeAtBBABBPM'` (`0.5E-3`)::
%	Threshold of change of trajectory offset at BBA-BPM, if below BBA
%	evaluation is not performed.
% `'minBPMrangeDownstream'` (`100E-6`)::
%	Minimum change of trajectory offset at downstream BPMs to be included
%	in calculation.
% `'maxKickatInjection'` (`0.9E-3`)::
%	 Maximum kick [rad] change at injection.
% `'maxOffsetatInjection'` (`0.9E-3`)::
%	 Maximum offset [m] change at injection.
% `'nSteps'` (`10`)::
%	 Number of different trajectories for each quadrupole setting
% `'nQvec'` (`[0.95 1.05]`)::
%	 Relative strength variation of quadrupoles.
% `'quadOrdPhaseAdvance'` (`[]`)::
%	 Quadrupole ordinate which is used to change the phase advance between
%	 injection and BBA-BPM.
% `'quadStrengthPhaseAdvance'` (`[0.95 0.8 1.05]`)::
%	 Relative quadrupole strength variation used to change the phase
%	 advance between injection and BBA-BPM.
% `'nXPointsNeededAtDownstreamBPM'` (`3`)::
%	 Number of x-positions at BBA-BPM required for linear regression at
%	 downstream BPMs.
% `'nXPointsNeededAtDownstreamBPM'` (`3`)::
%	 Number of x-positions at BBA-BPM required for linear regression at
%	 downstream BPMs.
% `'maxNumOfDownstreamBPMs'` (`length(SC.ORD.BPM)`)::
%	 Number downstream BPMs considered in the data evaluation.
% `'minSlopeForFit'` (`0.03`)::
%	 Minimum fitted slope at downstream BPMs with respect to quadrupole
%	 change which is still used in determining BPM center.
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
% *SCsetMags2SetPoints*, *SCgetBPMreading*



	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Input check

	p = inputParser;
	addOptional(p,'minBPMrangeAtBBABBPMTarget',1E-3);
	addOptional(p,'minBPMrangeAtBBABBPM',500E-6);
	addOptional(p,'minBPMrangeDownstream',100E-6);
	addOptional(p,'maxKickatInjection',.9E-3);
	addOptional(p,'maxOffsetatInjection',.9E-3);
	addOptional(p,'nSteps',10);
	addOptional(p,'nQvec',[0.95,1.05]);
	addOptional(p,'quadOrdPhaseAdvance',[ ]);
	addOptional(p,'quadStrengthPhaseAdvance',[0.95 0.8 1.05]);
	addOptional(p,'nXPointsNeededAtDownstreamBPM',3);
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

	% Allocate BBA results
	BPMcenter = nan(2,length(SC.ORD.BPM));
	BPMstd    = nan(2,length(SC.ORD.BPM));
	
	% Offset variation at injection
	kickVec0(1,:)  = par.maxOffsetatInjection * linspace(-1,1,par.nSteps);
	% Kick variation at injection
	kickVec0(2,:)  = par.maxKickatInjection   * linspace(-1,1,par.nSteps);

	% Save initial injection trajectory
	initialZ0 = SC.INJ.Z0;

	% Read quadrupole initial setpoint for changing the phase advance
	if ~isempty(par.quadOrdPhaseAdvance)
		quadPhaseAdvanceSP0 = SC.RING{par.quadOrdPhaseAdvance}.SetPointB(2);
	end
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Perform BBA routine

	% Loop over BPMs
	for jBPM=1:size(BPMords,2) % jBPM: Index of BPM adjacent to quad for BBA


		% Horizontal/vertical
		for nDim=1:2
			if par.verbose;fprintf('BBA-BPM %d/%d, nDim = %d',jBPM,size(BPMords,2),nDim);end

			
			% Get BPM index-ordinate pairing
			BPMind = find(BPMords(nDim,jBPM)==SC.ORD.BPM);
			
			% Define ordinate of quadrupole for BBA
			qOrd = QuadOrds(nDim,jBPM);
			% Read quadrupole initial setpoint
			qSP0 = SC.RING{qOrd}.SetPointB(2);


			% Scale maximum initial trajectory variation before beam gets lost
			[kickVec, BPMrange] = scaleInjectionToReachBPM(SC,BPMind,nDim,initialZ0,kickVec0,par);


			% Check if offset range at BBA-BPM is to small and vary phase advance
			if ~isempty(par.quadOrdPhaseAdvance) && BPMrange < par.minBPMrangeAtBBABBPMTarget
				[SC,kickVec] = scanPhaseAdvance(SC,BPMind,nDim,initialZ0,kickVec0,par);
			end

			% Perform data measurement
			[BPMpos,tmpTra] = dataMeasurement(SC,qOrd,BPMind,nDim,initialZ0,kickVec,par);

			% Reset BBA quad to nominal setpoint
			SC = SCsetMags2SetPoints(SC,qOrd,2,2,qSP0,...
				'method','abs',...
				'dipCompensation',1);
			% Reset phase advance quadrupole to initial setpoint
			if ~isempty(par.quadOrdPhaseAdvance)
				SC = SCsetMags2SetPoints(SC,par.quadOrdPhaseAdvance,2,2,quadPhaseAdvanceSP0,...
					'method','abs',...
					'dipCompensation',1);
			end

			% Perform data evaluation
			[BPMcenter(nDim,BPMind),BPMstd(nDim,BPMind)] = dataEvaluation(SC,BPMords,jBPM,BPMpos,tmpTra,nDim,qOrd,par);

			% Check if BPM center could be identified in measurement
			if ~isnan(BPMcenter(nDim,BPMind))
				% Apply BBA results
				SC.RING{BPMords(nDim,jBPM)}.Offset(nDim) = SC.RING{BPMords(nDim,jBPM)}.Offset(nDim) + BPMcenter(nDim,BPMind);
			end
		end

		% Plot BBA results
		if par.plotResults
			plotBBAResults(SC,jBPM,BPMords,BPMcenter,BPMstd)
		end
	end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data measurement
function [BPMpos,tmpTra] = dataMeasurement(SC,qOrd,BPMind,nDim,initialZ0,kickVec,par)

	% Prealocate x and y data
	tmpTra = nan(size(kickVec,2),length(par.nQvec),par.maxNumOfDownstreamBPMs);
	BPMpos = nan(size(kickVec,2),length(par.nQvec));

	% Loop over different quadrupole settings
	for nQ=1:length(par.nQvec)

		% Set quadrupole to different setpoints
		SC = SCsetMags2SetPoints(SC,qOrd,2,2,par.nQvec(nQ),...
			'method','rel',...
			'dipCompensation',1);

		% Loop over different trajectories
		for nKick=1:size(kickVec,2)

			% Vary initial trajectory
			SC.INJ.Z0(2*nDim)   = initialZ0(2*nDim  ) + kickVec(2,nKick); % kick angle
			SC.INJ.Z0(2*nDim-1) = initialZ0(2*nDim-1) + kickVec(1,nKick); % offset

			% Calculate beam reading
			B = SCgetBPMreading(SC);

			% JUST FOR PLOTTING
			if par.plotLines; plotBBAstep(SC,BPMind,nDim,nQ,qOrd,nKick,par);end

			% Save BBA-BPM readings
			BPMpos(nKick,nQ) = B(nDim,BPMind);

			% Save downstream BPM readings
			tmpTra(nKick,nQ,:) = B(nDim, (BPMind+1):(BPMind+par.maxNumOfDownstreamBPMs) );

			% Reset injected beam trajectory
			SC.INJ.Z0 = initialZ0;
		end
	end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot BBA step
function plotBBAstep(SC,BPMind,nDim,nQ,qOrd,nKick,par)
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
	figure(99);subplot(length(par.nQvec),1,nQ);hold on;
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data evaluation
function [BPMcenter,BPMstd] = dataEvaluation(SC,BPMords,jBPM,BPMpos,tmpTra,nDim,qOrd,par)
	if par.plotLines
% 		figure(56),clf;hold on; plot3(0,1E6*SC.RING{qOrd}.MagnetOffset(nDim),0,'rD','MarkerSize',40,'MarkerFaceColor','b'); hold on;
		figure(56),clf;hold on; plot3(0,1E6*SC.RING{qOrd}.T2(2*nDim-1),0,'rD','MarkerSize',40,'MarkerFaceColor','b'); hold on;
		set(gca,'ZLim',[-1 1],'XLim',[0 480],'box','on');
	end
 
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
			if length(x)>=par.nXPointsNeededAtDownstreamBPM ...        % More than 3 different points at BBA-BPM
					&& tmpRangeX(i) > par.minBPMrangeAtBBABBPM ...     % Range at BBA-BPM        > par.minBPMrangeAtBBABBPM
					&& tmpRangeY(i) > par.minBPMrangeDownstream        % Range at downstream BPM > pra.minBPMrangeDownstream

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
			if par.plotLines; plot3(repmat(jBPM+nBPM,size(x)),1E6*x,1E3*y,'ko');
				plot3(repmat(jBPM+nBPM,size(BPMpos)),1E6*BPMpos,1E3*(slope*BPMpos+offset),'k-')
				plot3(jBPM+nBPM,1E6*tmpCenter(nBPM),0,'Or','MarkerSize',10);
			end

		end
	end

	% Check if measurment requirements have been met
	if (max(tmpRangeX) < par.minBPMrangeAtBBABBPM) ...                % Max. range at BBA-BPM to small
			|| (max(tmpRangeY) < par.minBPMrangeDownstream) ...       % Max. range at downstream BPM to small
			|| std(tmpCenter,'omitnan') > par.maxStdForFittedCenters  % Fitted magnetic centers to far spread out

		BPMcenter = nan;
		BPMstd    = nan;
	else
		% Calculate mean magnetic center and error
		BPMcenter = mean(tmpCenter,'omitnan');
		BPMstd    = std( tmpCenter,'omitnan');
	end

	% Just plotting
	if par.plotLines; plot3(0,1E6*BPMcenter,0,'kD','MarkerSize',30,'MarkerFaceColor','r'); title(sprintf('BBA-BPM: %d // qOrd: %d // qFam: %s // nDim: %d // FinOffset = %.0f $\\mu m$',jBPM,qOrd,SC.RING{qOrd}.FamName,nDim,1E6*abs(SC.RING{BPMords(nDim,jBPM)}.Offset(nDim) + BPMcenter -SC.RING{qOrd}.MagnetOffset(nDim))));end
	if par.plotLines; xlabel('Index of BPM'); ylabel('BBA-BPM offset [$\mu$m]'); zlabel('Downstream offset [mm]');end
	if par.plotLines; set(findall(figure(56),'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');set(findall(gcf,'-property','FontSize'),'FontSize',28);set(gcf,'color','w');drawnow;end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale injection trajectory variation
function [kickVec, BPMrange] = scaleInjectionToReachBPM(SC,BPMind,nDim,initialZ0,kickVec0,par)
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
		if ~( BPMrange < par.minBPMrangeAtBBABBPMTarget )
			break
		end
	end

	% Check if offset range at BBA-BPM is still to small
	if BPMrange < par.minBPMrangeAtBBABBPMTarget
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
% Plot BBA results
function plotBBAResults(SC,jBPM,BPMords,BPMcenter,BPMstd)
	figure(90),clf;subplot(2,3,1:3)
	% Define figure of merit
	tmpBBAchange = BPMcenter;
	tmpBBAchange(isnan(tmpBBAchange))=0;
	fom  = cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'Offset'))';
	fom0 = cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'Offset'))' - tmpBBAchange;
	fom(isnan(BPMstd))=nan;

	nHist = 1.1*max(abs(fom0(:)))*linspace(-1,1,floor(length(SC.ORD.BPM)/2));
	

	% Plot histogram of new BPM offsets
	[offset,slope] = hist(fom',nHist);
	plot(1E6*repmat(slope,1,size(offset,2)),offset,'LineWidth',2);hold on
	xlabel('Remaining BPM offset [$\mu$m]'),ylabel('Number of counts');
	% Plot histogram of initial BPM offsets
	[offset,slope] = hist(fom0(:),nHist);
	plot(1E6*slope,offset,'k-','LineWidth',2)
	legend({sprintf('Horizontal rms: $%.0f\\mu m$',1E6*sqrt(mean(fom(1,:).^2,'omitnan'))),sprintf('Vertical rms: $%.0f\\mu m$',1E6*sqrt(mean(fom(2,:).^2,'omitnan'))),sprintf('Uncorrected rms: $%.0f\\mu m$',1E6*sqrt(mean(fom0(:).^2,'omitnan')))},'Interpreter','latex')

	subplot(2,3,4:5)
	% Plot error bar of new BPM offset
	errorbar(1E6*abs(fom'),1E6*abs(BPMstd'),'O','LineWidth',2)
	ylabel('Remaining BPM offset [$\mu$m]'),xlabel('Index of BPM');set(gca,'XLim',[1 length(SC.ORD.BPM)])
	legend({sprintf('Successfull BPMs: %d/%d',length(find(~isnan(abs(fom(1,:))))),length(BPMords(1,1:jBPM))),sprintf('Successfull BPMs: %d/%d',length(find(~isnan(abs(fom(2,:))))),length(BPMords(1,1:jBPM)))})

	subplot(2,3,6)
	% Plot new BPM offset vs estimated BPM offset fit error
	plot(1E6*abs(fom'),1E6*abs(BPMstd'),'O','LineWidth',2);hold on
	plot([0 1E6*max(max([abs(fom(:));abs(BPMstd(:))]))],[0 1E6*max(max([abs(fom(:));abs(BPMstd(:))]))],'k--')
	xlabel('Remaining BPM offset [$\mu$m]'),ylabel('Offset error [$\mu$m]');

	set(gcf,'color','w');set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');set(findall(gcf,'-property','FontSize'),'FontSize',18);
	drawnow
end
