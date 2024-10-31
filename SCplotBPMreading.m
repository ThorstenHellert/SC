function SCplotBPMreading(SC,B,T)
% SCplotBPMreading
% ================
%
% NAME
% ----
% SCplotBPMreading - plots BPM readings and particle trajectories
%
% SYNOPSIS
% --------
% `SCplotBPMreading(SC, B, T)`
%
%
% DESCRIPTION
% -----------
% `SCplotBPMreading` is a plotting function used by *SCgetBPMreading* and plots the BPM readings and
% particle trajectories including the aperture model. Note that the BPM readings must include all
% registered BPMs and the trajectories must be evaluated at all lattice elements.
% 
%
% INPUTS
% ------
% `SC`:: SC base structure
% `B`::  BPM readings
% `T`::  Trajectories
%
% SEE ALSO
% --------
% *SCgetBPMreading*


	ApertureForPLotting = getRingAperture(SC);

	% Check if orbit mode
	if strcmp(SC.INJ.trackMode,'ORB')
		SC.INJ.nTurns     = 1;
		SC.INJ.nParticles = 1;
	end
	figure(23);clf;	tmpCol = get(gca, 'ColorOrder');
	% Labels and legend
	ylabelStr = {'$\Delta x$ [mm]','$\Delta y$ [mm]'};
	legStr = {'Particle trajectories','BPM reading','Aperture'};

	% s positions of all lattice elements
	sPos = findspos(SC.RING,1:length(SC.RING))';
	sMax = findspos(SC.RING,length(SC.RING)+1);

	% Loop over dimensions
	for nDim=1:2
		ax(nDim)=subplot(2,1,nDim);hold on

		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Trajectories
		x = reshape([0:SC.INJ.nTurns-1]*sMax+sPos,1,[]);
		x = repmat(x',1,SC.INJ.nParticles);
		for nS=1:SC.INJ.nShots
			M = SCparticlesIn3D(T(:,:,nS),SC.INJ.nParticles);
			y = 1E3*squeeze(M(2*nDim-1,:,:));
			legVec=plot(x,y,'k');
		end
		% Use one line for legend handle
		legVec=legVec(1);

		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% BPMs
		x    = reshape([0:SC.INJ.nTurns-1]*sMax+findspos(SC.RING,SC.ORD.BPM)',1,[]);
		y    = 1E3*(B(nDim,:));
		legVec(2) = plot(x,y,'rO');

		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Aperture
		if ~isempty(ApertureForPLotting)
			apS = sPos(ApertureForPLotting.apOrds);

			x = reshape([0:SC.INJ.nTurns-1]*sMax+apS,1,[]);
			y = 1E3*repmat(ApertureForPLotting.apVals{nDim},1,SC.INJ.nTurns);

			legVec(3)=stairs(x,y(1,:),'Color',tmpCol(1,:),'LineWidth',4);
			stairs(x,y(2,:),'Color',tmpCol(1,:),'LineWidth',4)
		end

		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Ring ending
		x = reshape([0:SC.INJ.nTurns-1]*sMax,1,[]);
		for nT=1:SC.INJ.nTurns
			plot(x(nT)*[1 1],10*[-1 1],'k:')
		end

		% Set axes properties
		set(gca,...
			'XLim',[0 SC.INJ.nTurns*sPos(end)],...
			'box','on',...
			'Position',get(gca,'Position')+[0 .07*(nDim-1) 0 0])

		set(gca,'ylim',.5*[-1 1])
		
		% Label and legend
		ylabel(ylabelStr{nDim});
		legend(legVec,legStr{1:length(legVec)});

	end

	% Label, white background, font size and latex
	xlabel('$s$ [m]');
	set(findall(gcf,'-property','FontSize'),'FontSize',18);
	set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
	set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
	set(gcf,'color','w');
	drawnow
		
	% Link x-axis
	linkaxes([ax(1) ax(2)],'x')
	

	
%   % (needed for making movies)
% 	global figFrame 
% 	tmp = getframe(gcf);
% 	figFrame(end+1).cdata    = tmp.cdata;
% 	figFrame(end).colormap = tmp.colormap;

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions

% Get aperture
function ApertureForPLotting = getRingAperture(SC)
	apOrds=[];apVals={[],[]};
	for ord=1:length(SC.RING)
		if isfield(SC.RING{ord},'EApertures') || isfield(SC.RING{ord},'RApertures')
			apOrds(end+1) = ord;
			for nDim=1:2
				if isfield(SC.RING{ord},'EApertures')
					apVals{nDim}(:,end+1) = SC.RING{ord}.EApertures(nDim) * [-1 1];
				else
					apVals{nDim}(:,end+1) = SC.RING{ord}.RApertures(2*(nDim-1) + [1 2] );
				end
			end
		end
	end
	if isempty(apOrds)
		ApertureForPLotting = [];
	else
		ApertureForPLotting.apOrds = apOrds;
		ApertureForPLotting.apVals = apVals;
	end
end
