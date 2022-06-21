function SCplotPhaseSpace(SC,varargin)
% SCplotPhaseSpace
% ================
%
% NAME
% ----
% SCplotPhaseSpace - Plots the turn-by-turn phase space
%
% SYNOPSIS
% --------
% `SCplotPhaseSpace(SC [,options])`
%
%
% DESCRIPTION
% -----------
% Plots the horizontal, vertical and longitudinal turn-by-turn phase space for a bunch of partciles
% as specified in the injection pattern `SC.INJ`. If `'RFCavityPass'` is the specified cavity pass
% method, the slippage length per turn caused by a potential frequency error is calculated and the
% longitudinal coordinate is adjusted for plotting.
%
%
% INPUTS
% ------
% `SC`:: SC base structure
%
%
% OPTIONS
% -------
% The following options can be given as name-vale pairs:
%
% `'nTurns'` (`SC.INJ.nTurns`)::
%   Number of turns.
% `'nParticles'` (`SC.INJ.nParticles`)::
%   Number of particles.
% `'ord'` (1)::
%   Ordinate for plotting.
% `'plotCO'` (0)::
%   Logical value defining if closed orbit should be plotted as well.
% `'customBunch'` ([])::
%   [6 x N] array of particle start points for phase space evaluation. If not specified a bunch is
%   generated using *SCgenBunches*.


% Parse input
p = inputParser;
addOptional(p,'ord',1);
addOptional(p,'plotCO',0);
addOptional(p,'costumBunch',[]); % Typo in old version, optional argument will be removed soon
addOptional(p,'customBunch',[]);
addOptional(p,'nParticles',SC.INJ.nParticles);
addOptional(p,'nTurns',SC.INJ.nTurns);
parse(p,varargin{:});
par = p.Results;

if length(par.ord)~=1
	error('Not implemented anymore...')
end

% Adjust number of particles/turns
SC.INJ.nParticles = par.nParticles;
SC.INJ.nTurns     = par.nTurns;

% Deal with typo in old version
if isempty(par.customBunch) && ~isempty(par.costumBunch)
	warning('Typo in old version, please use ''customBunch'' going forward')
	par.customBunch = par.costumBunch;
end

% Generate initial particle distribution
if isempty(par.customBunch)
	Zin = SCgenBunches(SC);
else
	Zin = par.customBunch;
	SC.INJ.nParticles = size(Zin,2);
end

% Calculate trajectory
T = atpass(SC.RING, Zin, 1, SC.INJ.nTurns, par.ord);
T(:,isnan(T(1,:))) = nan;
T3D = SCparticlesIn3D(T,SC.INJ.nParticles);

% Define labels and title
labelStr={'$\Delta x$ [$\mu$m]','$\Delta x''$ [$\mu$rad]','$\Delta y$ [$\mu$m]','$\Delta y''$ [$\mu$rad]','$\Delta S$ [m]','$\delta E$ $[\%]$'};
titleStr = {'Horizontal','Vertical','Longitudinal'};

% Adjust phase coordinate so it not runs away
if strcmp(SC.RING{SC.ORD.Cavity(1)}.PassMethod,'RFCavityPass')
	L0_tot=0;
	for i=1:length(SC.RING)
		L0_tot=L0_tot+SC.RING{i}.Length;
	end
	% Calculate slippage as a function of turns
	lengthSlippage = 299792458 * (SC.RING{SC.ORD.Cavity(1)}.HarmNumber/SC.RING{SC.ORD.Cavity(1)}.Frequency - L0_tot/299792458);
	% Adjust s coordinate
	T3D(6,:,:) = T3D(6,:,:) - lengthSlippage * [1:par.nTurns];

	labelStr{5} = '$\Delta S_{act}$ [m]';
end

% Calculate closed orbit
if par.plotCO
	CO = findorbit6(SC.RING,par.ord);
	if isnan(CO(1))
		startPointGuess = mean(squeeze(mean(T3D,2,'omitnan')),2,'omitnan');
		CO = findorbit6(SC.RING,par.ord,startPointGuess);
		if isnan(CO(1))
			CO = nan(6,1);
		end
	end
else
	CO = nan(6,1);
end

% Normalize to plotting units
T3D       = T3D .* repmat([1E6;1E6;1E6;1E6;1E2;1],1,par.nTurns);
SC.INJ.Z0 = SC.INJ.Z0 .* [1E6;1E6;1E6;1E6;1E2;1];
CO        = CO  .* repmat([1E6;1E6;1E6;1E6;1E2;1],1:length(par.ord));

% Switching 5th and 6th coordinate
T3D([5 6],:)     = T3D([6 5],:);
CO([5 6],:)      = CO([6 5],:);
SC.INJ.Z0([5 6]) = SC.INJ.Z0([6 5]);


% Create figure
figure(100);clf;pVec=[];legStr=[];
for nType=1:3
	subplot(1,3,nType);hold on

	% Particles
	for nP=1:SC.INJ.nParticles
		x = T3D(2*nType-1,:,nP);
		y = T3D(2*nType  ,:,nP);
		scatter(x,y,10,1:SC.INJ.nTurns)
	end

	% Injected beam centroid
	pVec(1)=plot(SC.INJ.Z0(2*nType-1), SC.INJ.Z0(2*nType),'O','MarkerSize',15,'LineWidth',3);
	legStr{1} = 'Injection point';

	% Closed orbit
	if par.plotCO
		pVec(2)=plot(CO(2*nType-1), CO(2*nType),'X','MarkerSize',20,'LineWidth',3);
		legStr{2} = 'Closed orbit';
	end

	% Labels and title
	set(gca,'box','on');
	xlabel(labelStr{2*nType-1});
	ylabel(labelStr{2*nType},'Interpreter','none');
	title(sprintf('%s @Ord: %d',titleStr{nType},par.ord))
end

% Legend
legend(pVec,legStr)

% Colorbar to show turns
c=colorbar;
c.Label.String = 'Number of turns';
c.Label.Interpreter = 'Latex';


% Make chic
set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(gcf,'color','w');drawnow

end
