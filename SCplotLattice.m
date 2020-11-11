function SCplotLattice(SC,varargin)
% SCplotLattice
% =============
%
% NAME
% ----
% SCplotLattice - Plots the lattice including the location of BPMs and CMs
%
% SYNOPSIS
% --------
% `SCplotLattice(SC [, options])`
%
%
% DESCRIPTION
% -----------
% This function plots the lattice functions as well as the distribution of registered magnets 
% including CMs and BPMs.
%
%
% INPUTS
% ------
% `SC`:: SC base structure
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'transferLine'` (0)::
%   If true the function 'twissline' is used to calculate the lattice functions
% `'sRange'` ([])::
%   Array ['sMin','sMax'] defining the plot range [m].
% `'oList'` ([])::
%   If `'sRange'` is empty, `'oList'`  can be used to specify a list of ordinates at which the lattice should be plotted
% `'nSectors'` (1)::
%   If `'oList'` is empty, `'nSectors'` can be used to plot only the first fraction of the lattice
% `'plotIdealRing'` (1)::
%   Specify if 'SC.IDEALRING' should be used to plot twiss functions, otherwise 'SC.RING'.
% `'fontSize'` (16)::
%   Figure font size.
%
% EXAMPLES
% --------
% Plots the complete lattice for a ring
% ------------------------------------------------------------------
% SCplotLattice(SC);
% ------------------------------------------------------------------
%
% Plots the lattice from ordinate 30 to 130 for a transfer line
% ------------------------------------------------------------------
% SCplotLattice(SC,'transferLine',1,'oList',30:120);
% ------------------------------------------------------------------
% 
% Plots the lattice of one arc for a twelfe-fold symmetric ring lattice
% ------------------------------------------------------------------
% SCplotLattice(SC,'nSectors',12);
% ------------------------------------------------------------------

	% Parse input
	p = inputParser;
	p.KeepUnmatched=1;
	addOptional(p,'transferLine',0);
	addOptional(p,'nSectors',1);
	addOptional(p,'oList',[]);
	addOptional(p,'plotIdealRing',1);
	addOptional(p,'sRange',[]);
	addOptional(p,'fontSize',16);
	parse(p,varargin{:});
	par = p.Results;


	% Get s-positions along the lattice
	sPos = findspos(SC.RING,1:length(SC.RING))';

	% Check if oList is given explicitly
	if isempty(par.oList)
		% Take lattice elements smaller than 1/nSectors of total length
		par.oList = find(sPos<=(sPos(end)/par.nSectors))';
	end
	
	% Check if range of s-positions is given
	if ~isempty(par.sRange)
		par.oList = intersect(find(sPos>=par.sRange(1)),find(sPos<=par.sRange(2)))';
	end
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Check for transfer line or for ring
	if par.transferLine
		% Check if initial Courant-Snyder paramaters are defined in transfer line lattice
		if ~isfield(SC.RING{1},'TD')
			warning('Transfer line lattice did not contain initial parameters needed for beta function calculation!')
			beta = nan(2,length(par.oList));
			disp = nan(1,length(par.oList));
		else
			% Get beta functions and dispersion
			TD   = twissline(SC.IDEALRING,0,SC.IDEALRING{1}.TD,par.oList,'chrom',1E-8);
			beta = reshape([TD.beta],2,[]);
			disp = reshape([TD.Dispersion],4,[]);
		end
	else
		% Get beta functions and dispersion
		if par.plotIdealRing
			[ld,~,~] = atlinopt(SC.IDEALRING,1e-3,par.oList);
		else
			[ld,~,~] = atlinopt(SC.RING,1e-3,par.oList);
		end

		% Generate arrays from arkward linopt output
		beta = reshape([ld.beta],2,[]);
		disp = reshape([ld.Dispersion],4,[]);
	end

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Find magnets in lattice
	DIP=[];QUAD=[];SEXT=[];OCT=[];SKEW=[];
	for ord=par.oList
		% Search for dipoles
		if isfield(SC.RING{ord},'BendingAngle') && SC.RING{ord}.BendingAngle~=0
			DIP(end+1) = ord;
		end
		if isfield(SC.RING{ord},'NomPolynomB')
			% Search for nominal quadrupole fields
			if any(find(SC.RING{ord}.NomPolynomB)==2)
				QUAD(end+1) = ord;
			end
			% Search for nominal sextupole fields
			if any(find(SC.RING{ord}.NomPolynomB)==3)
				SEXT(end+1) = ord;
			end
			% Search for nominal octupole fields
			if any(find(SC.RING{ord}.NomPolynomB)==4)
				OCT(end+1) = ord;
			end
		end
	end
	if isfield(SC,'ORD') && isfield(SC.ORD,'SkewQuad')
		SKEW = SC.ORD.SkewQuad(ismember(SC.ORD.SkewQuad,par.oList));
	end

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get aperture values from lattice structure
	ApertureForPlotting.apOrds=[];ApertureForPlotting.apVals={[],[]};
	for nEl=par.oList
		if isfield(SC.RING{nEl},'EApertures') || isfield(SC.RING{nEl},'RApertures')
			ApertureForPlotting.apOrds(end+1) = nEl;
			for nDim=1:2
				if isfield(SC.RING{nEl},'EApertures')
					ApertureForPlotting.apVals{nDim}(:,end+1) = SC.RING{nEl}.EApertures(nDim)*[-1 1];
				else
					ApertureForPlotting.apVals{nDim}(:,end+1) = SC.RING{nEl}.RApertures( 2*(nDim-1) + [1 2] );
				end
			end
		end
	end

% 	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	% Get support structures (usually to messy if included, see 2nd subplot)
% 	GIRDER=[];PLINTH=[];SECTION=[];
%
% 	if isfield(SC.ORD,'Girder')
% 		GIRDER = SC.ORD.Girder(:,2==(SC.ORD.Girder(1,:)>=par.oList(1)) + (SC.ORD.Girder(2,:)<=par.oList(end)));
% 	end
% 	if isfield(SC.ORD,'Plinth')
% 		PLINTH = SC.ORD.Plinth(:,2==(SC.ORD.Plinth(1,:)>=par.oList(1)) + (SC.ORD.Plinth(2,:)<=par.oList(end)));
% 	end
% 	if isfield(SC.ORD,'Section')
% 		SECTION = SC.ORD.Section(:,2==(SC.ORD.Section(1,:)>=par.oList(1)) + (SC.ORD.Section(2,:)<=par.oList(end)));
% 	end

	% Create figure
	figure(1);clf;hold on;tmpCol=get(gca, 'ColorOrder');



	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot twiss functions
	ax1=subplot(3,1,1);hold on;
	% PLot dispersion
	yyaxis right
	plot(sPos(par.oList),1E2*disp(1,:),'LineWidth',2)
	ylabel('$\eta_x$ [cm]');
	% Plot beta functions
	yyaxis left
	plot(sPos(par.oList),beta,'LineWidth',2)
	ylabel('$\beta$ [m]');

	% Create some space and generate the legend
	set(gca,'ylim',get(gca,'ylim').*[1 1.5],'xlim',[sPos(min(par.oList)) sPos(max(par.oList))],'box','on')
	legend({'Hor. Beta','Ver. Beta','Hor. Disp.'},'Location','N','Orientation','Horizontal');title('Beta Functions and Dispersion');%xlabel('$s$ [m]');


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot quadrupoles, dipoles and sextupoles
	ax2=subplot(3,1,2);hold on;legStr={};legVec=[];
	title('Aperture and Magnets');ylabel('Aperture [mm]');set(gca,'xlim',[sPos(min(par.oList)) sPos(max(par.oList))],'box','on');%xlabel('$s$ [m]')

	% Aperture % %% %% %% %
	if ~isempty(ApertureForPlotting.apOrds)
		apS = sPos(ApertureForPlotting.apOrds);
		lStyle = {'-',':'};
		for nDim=1:2
			p3(nDim)=stairs(apS,1E3*ApertureForPlotting.apVals{nDim}(1,:),'Color',tmpCol(nDim,:),'LineWidth',4,'LineStyle',lStyle{nDim});
			stairs(apS,1E3*ApertureForPlotting.apVals{nDim}(2,:),'Color',tmpCol(nDim,:),'LineWidth',4,'LineStyle',lStyle{nDim});
		end
		legVec=p3;legStr={'Hor. Ap.','Ver. Ap.'};
	end
	
	% Get scaling factor for magent size
	scale = 1E3*max([max(abs(ApertureForPlotting.apVals{1}(:))) max(abs(ApertureForPlotting.apVals{2}(:)))]);
	if isempty(scale)
		scale = 10;
	end
	
	% Draw octupole magnets
	for nM=1:length(OCT)
		rectangle('Position',[sPos(OCT(nM)),-scale,sPos(OCT(nM)+1)-sPos(OCT(nM)),scale ],'FaceColor',tmpCol(6,:));
		if nM==1;legStr{end+1}='Oct';legVec(end+1)=bar(-1,0,'FaceColor',tmpCol(6,:));end
	end
	% Draw sextupole magnets
	for nM=1:length(SEXT)
		rectangle('Position',[sPos(SEXT(nM)),0,sPos(SEXT(nM)+1)-sPos(SEXT(nM)),scale ],'FaceColor',tmpCol(5,:));
		if nM==1;legStr{end+1}='Sext';legVec(end+1)=bar(-1,0,'FaceColor',tmpCol(5,:));end
	end
	% Draw dipole magnets
	for nM=1:length(DIP)
		rectangle('Position',[sPos(DIP(nM)),0,sPos(DIP(nM)+1)-sPos(DIP(nM)),scale/2 ],'FaceColor',[0 0 0]);
		if nM==1;legStr{end+1}='Dip';legVec(end+1)=bar(-1,0,'k');end
	end
	% Draw quadrupole magnets
	for nM=1:length(QUAD)
		rectangle('Position',[sPos(QUAD(nM)),-scale/2,sPos(QUAD(nM)+1)-sPos(QUAD(nM)),scale/2 ],'FaceColor',tmpCol(3,:));
		if nM==1;legStr{end+1}='Quad';legVec(end+1)=bar(-1,0,'FaceColor',tmpCol(3,:));end
	end

% 	% Draw support structures (usually to messy if included)
% 	tmp = get(gca,'ylim');colG={'k','r'};colP={'k',tmpCol(6,:)};colS={'k',tmpCol(7,:)};
% 	for n=1:size(GIRDER,2)
% 		plot([sPos(GIRDER(1,n)) sPos(GIRDER(2,n))],1.1*tmp(1)*[1 1],'Color',colG{1+mod(n,2)},'LineWidth',2)
% 	end
% 	for n=1:size(PLINTH,2)
% 		plot([sPos(PLINTH(1,n)) sPos(PLINTH(2,n))],1.2*tmp(1)*[1 1],'Color',colP{1+mod(n,2)},'LineWidth',2)
% 	end
% 	for n=1:size(SECTION,2)
% 		plot([sPos(SECTION(1,n)) sPos(SECTION(2,n))],1.3*tmp(1)*[1 1],'Color',colS{1+mod(n,2)},'LineWidth',2,'LineStyle','--')
% 	end


	% Create some space and generate the legend
	set(gca,'ylim',get(gca,'ylim').*[1 1.3])
	legend(legVec,legStr,'Location','North','Orientation','Horizontal')


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot BPMs, CMs and skew quads
	ax3=subplot(3,1,3);hold on;legStr={};legVec=[];

	% Draw CMs
	if isfield(SC,'ORD') && isfield(SC.ORD,'CM')
		for nDim=1:2
			for ord=intersect(SC.ORD.CM{nDim},par.oList)
				if strcmp(SC.RING{ord}.PassMethod,'CorrectorPass')
					legVec(nDim) = bar(sPos(ord),(-1)^(nDim-1)*4,'FaceColor',tmpCol(nDim,:),'BarWidth',(max(sPos(par.oList))-min(sPos(par.oList)))/100);
				else
					rectangle('Position',[sPos(ord),0-4*(nDim-1),sPos(ord+1)-sPos(ord),4],'FaceColor',tmpCol(nDim,:));
					% Fake plot for legend
					legVec(nDim) = bar(-1,0,'FaceColor',tmpCol(nDim,:));
				end
			end
		end
		legStr = {'HCM','VCM'};
	end

	% Draw skew quadrupole magnets
	for nM=1:length(SKEW)
		rectangle('Position',[sPos(SKEW(nM)),-2,sPos(SKEW(nM)+1)-sPos(SKEW(nM)),4 ],'FaceColor',tmpCol(5,:));
		if nM==1;legStr{end+1}='SKEW';legVec(end+1)=bar(-1,0,'FaceColor',tmpCol(5,:));end
	end


	% Draw BPMs
	if isfield(SC,'ORD') && isfield(SC.ORD,'BPM')
		for ord=intersect(SC.ORD.BPM,par.oList)
			rectangle('Position',[sPos(ord)-diff(sPos(par.oList([1 end])))/300,-3,diff(sPos(par.oList([1 end])))/150,6],'FaceColor','k');
		end
		% Fake plot for legend
		legVec(end+1) = bar(-1,0,'FaceColor','k');legStr{end+1}='BPM';
	end

	% Create some space and generate the legend
	set(gca,'ylim',get(gca,'ylim').*[1 1.3])
	legend(legVec,legStr,'Location','N','Orientation','Horizontal')
	title('BPMs and CMs');xlabel('$s$ [m]');set(gca,'box','on','xlim',[sPos(min(par.oList)) sPos(max(par.oList))],'yTickLabel','')

	% Link x-axis
	linkaxes([ax1 ax2 ax3],'x')
	
	% Make nice
	set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
	set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
	set(findall(gcf,'-property','FontSize'),'FontSize',par.fontSize);
	set(gcf,'color','w');
	drawnow
	%  set(gcf,'Position',[2284 709 1573 529]);
	%  print('-clipboard','-dbitmap','-r300')


