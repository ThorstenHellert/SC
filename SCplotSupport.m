function SCplotSupport(SC,varargin)
% SCplotSupport
% =============
%
% NAME
% ----
% SCplotSupport - Plots the offset and rolls of magnets, the support structure and BPMs
%
% SYNOPSIS
% --------
% `SCplotSupport(SC)`
%
%
% DESCRIPTION
% -----------
% This function plots the overall offsets and rolls of all magnets and BPMs, as well as the individual contributions
% from different support structures (if registered).
%
%
% INPUTS
% ------
% `SC`:: SC base structure
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'fontSize'` (16)::
%   Figure font size.
%
%
% SEE ALSO
% --------
% *SCgetSupportOffset*, *SCupdateSupport*

	% Parse input
	p = inputParser;
	addOptional(p,'fontSize',16);
	parse(p,varargin{:});
	par = p.Results;

	% Check if magnets and BPMs are registered
	if ~isfield(SC.ORD,'Magnet')
		error('Magnets must be registered. Use ''SCregisterMagnets''.')
	elseif ~isfield(SC.ORD,'BPM')
		error('BPMs must be registered. Use ''SCregisterBPMs''.')
	end

	% Create figure
	figure(1213);clf;tmpCol=get(gca, 'ColorOrder');ax=[];
	yLabStr = {'$\Delta x$ [$\mu$m]','$\Delta y$ [$\mu$m]'};

	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Offsets

	
	% Loop over both transverse planes
	for nDim=1:2
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot BPM offset
		ax(end+1)=subplot(12,1,4*(nDim-1)+ 1);hold on

		sBPM         = findspos(SC.RING,SC.ORD.BPM);
		offBPM       = cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'Offset'));
		offBPMStruct = cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'SupportOffset'));
		
		
		% Plot random BPM offset
		plot(sBPM,1E6*offBPM(:,nDim),'O','Color',tmpCol(2,:),'MarkerSize',6);
		% Plot BPM offset from support structure
		plot(sBPM,1E6*offBPMStruct(:,nDim),'-','Color',tmpCol(2,:));

		% Legend and axis stuff
		legend({'Random BPM offset','BPM support offset'});
		set(gca,'xlim',[0 findspos(SC.RING,length(SC.RING)+1)],'box','on','XTickLabel','')
		ylabel(yLabStr{nDim});

		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot overall magnet offset
		ax(end+1)=subplot(12,1,4*(nDim-1)+ 2);

		offMag=[];
		% Loop over magnets
		for ord=SC.ORD.Magnet
			% Get overall magnet offset
			offMag(1:2,end+1)=SC.RING{ord}.T2([1 3]);
		end
		% Plot magnet offset
		plot(findspos(SC.RING,SC.ORD.Magnet),1E6*offMag(nDim,:),'kO-');

		% Legend and axis stuff
		legend('Overall magnet offset');
		set(gca,'xlim',[0 findspos(SC.RING,length(SC.RING)+1)],'box','on','XTickLabel','')
		ylabel(yLabStr{nDim});

		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot individual offset contributions
		ax(end+1)=subplot(12,1,4*(nDim-1)+ [3:4]);hold on;

		% Get s position and support structure offset
		C = findspos(SC.RING,length(SC.RING)+1);
		s = linspace(0,C,100*length(SC.RING));
		off = SCgetSupportOffset(SC,s);

		% Plot line of support offset
		stairs(s,1E6*off(nDim,:));
		pVec=[];legStr=[];

		offMagSupport=[];offMag=[];i=0;
		% Loop over magnets
		for ord=SC.ORD.Magnet
			i=i+1;
			% Check if support structure is registered
			if isfield(SC.RING{ord},'SupportOffset')
				offMagSupport(1:2,i)=SC.RING{ord}.SupportOffset;
			else
				offMagSupport(1:2,i)=nan(1,2);
			end
			% Get magnet-to-magnet offset
			offMag(1:2,i)=SC.RING{ord}.MagnetOffset;
		end
		% Plot support offset
		plot(findspos(SC.RING,SC.ORD.Magnet),1E6*offMagSupport(nDim,:),'D','Color',tmpCol(1,:));
		% Fake plot for legend
		pVec(end+1)=plot([-2 -1],[0 0],'-D','Color',tmpCol(1,:));legStr{end+1}='Overall support offset';
		% Plot magnet offset
		pVec(end+1)=plot(findspos(SC.RING,SC.ORD.Magnet),1E6*offMag(nDim,:),'kx','MarkerSize',8);legStr{end+1}='Magnet offset';



		% Check if plinths are registered
		if isfield(SC.ORD,'Plinth')
			offPa=[];offPb=[];
			% Get plinth start offsets
			for ord=SC.ORD.Plinth(1,:)
				offPa(1:2,end+1)=SC.RING{ord}.PlinthOffset(:);
			end
			% Get plinth ending offsets
			for ord=SC.ORD.Plinth(2,:)
				offPb(1:2,end+1)=SC.RING{ord}.PlinthOffset(:);
			end

			% Plot plinth offset
			for i=1:size(SC.ORD.Plinth,2)
				plot(findspos(SC.RING,SC.ORD.Plinth(:,i)),1E6*[offPa(nDim,i),offPb(nDim,i)],'Color','k','LineWidth',4)
			end
			% Fake plot for legend entry
			pVec(end+1)=plot([-2 -1],[0 0],'Color','k','LineWidth',4);legStr{end+1}='Plinth offset';
		end
		% Check if girders are registered
		if isfield(SC.ORD,'Girder')
			offGa=[];offGb=[];
			% Get girder start offsets
			for ord=SC.ORD.Girder(1,:)
				offGa(1:2,end+1)=SC.RING{ord}.GirderOffset(:);
			end
			% Get girder ending offsets
			for ord=SC.ORD.Girder(2,:)
				offGb(1:2,end+1)=SC.RING{ord}.GirderOffset(:);
			end

			% Plot girder offset
			for i=1:size(SC.ORD.Girder,2)
				plot(findspos(SC.RING,SC.ORD.Girder(:,i)),1E6*[offGa(nDim,i),offGb(nDim,i)],'Color',tmpCol(5,:),'LineWidth',4)
			end
			% Fake plot for legend entry
			pVec(end+1)=plot([-2 -1],[0 0],'Color',tmpCol(5,:),'LineWidth',4);legStr{end+1}='Girder offset';
		end
		% Check if sections are registered
		if isfield(SC.ORD,'Section')
			offSa=[];offSb=[];
			% Get section start offsets
			for ord=SC.ORD.Section(1,:)
				offSa(1:2,end+1)=SC.RING{ord}.SectionOffset(:);
			end
			% Get section ending offsets
			for ord=SC.ORD.Section(2,:)
				offSb(1:2,end+1)=SC.RING{ord}.SectionOffset(:);
			end

			% Plot section offset
			for i=1:size(SC.ORD.Section,2)
				plot(findspos(SC.RING,SC.ORD.Section(:,i)),1E6*[offSa(nDim,i),offSb(nDim,i)],'Color',tmpCol(7,:),'LineWidth',3,'LineStyle',':')
			end
			% Fake plot for legend entry
			pVec(end+1)=plot([-2 -1],[0 0],'Color',tmpCol(7,:),'LineWidth',3,'LineStyle',':');legStr{end+1}='Section offset';
		end

		% Legend and axis stuff
		legend(pVec,legStr);
		set(gca,'xlim',[0 findspos(SC.RING,length(SC.RING)+1)],'box','on','XTickLabel','')
		ylabel(yLabStr{nDim});
	end





	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Rolls
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



	yLabStr = {'$\Delta \Phi$ [$\mu$rad]'};
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot BPM Roll
	ax(end+1)=subplot(12,1,9);hold on

	sBPM          = findspos(SC.RING,SC.ORD.BPM);
	rollBPM       = atgetfieldvalues(SC.RING(SC.ORD.BPM),'Roll');
	rollBPMStruct = atgetfieldvalues(SC.RING(SC.ORD.BPM),'SupportRoll');

	% Plot random BPM roll
	plot(sBPM,1E6*rollBPM,'O','Color',tmpCol(2,:),'MarkerSize',6);
	% Plot BPM roll from support structure
	plot(sBPM,1E6*rollBPMStruct,'-','Color',tmpCol(2,:));

	% Legend and axis stuff
	legend({'Random BPM roll','BPM support roll'});
	set(gca,'xlim',[0 findspos(SC.RING,length(SC.RING)+1)],'box','on','XTickLabel','')
	ylabel(yLabStr);

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot overall magnet roll
	ax(end+1)=subplot(12,1,10);

	rollMag=[];
	% Loop over magnets
	for ord=SC.ORD.Magnet
		% Get overall magnet roll
		rollMag(end+1)=SC.RING{ord}.RollAngle;
	end
	% Plot magnet roll
	plot(findspos(SC.RING,SC.ORD.Magnet),1E6*rollMag,'kO-');

	% Legend and axis stuff
	legend('Overall magnet roll');
	set(gca,'xlim',[0 findspos(SC.RING,length(SC.RING)+1)],'box','on','XTickLabel','')
	ylabel(yLabStr);


	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot individual roll contributions
	ax(end+1)=subplot(12,1,[11 12]);hold on;

	% Get s position and support structure roll
	C = findspos(SC.RING,length(SC.RING)+1);
	s = linspace(0,C,100*length(SC.RING));
	off = SCgetSupportOffset(SC,s);

	% Plot line of support roll
	pVec=[];legStr=[];

	rollMagSupport=[];rollMag=[];i=1;
	% Loop over magnets
	for ord=SC.ORD.Magnet
		% Check if support structure is registered
		if isfield(SC.RING{ord},'SupportRoll')
			rollMagSupport(i)=SC.RING{ord}.SupportRoll;
		else
			rollMagSupport(i)=nan;
		end
		% Get magnet-to-magnet roll
		rollMag(i)=SC.RING{ord}.MagnetRoll;
		i=i+1;
	end
	% Plot support roll
	plot(findspos(SC.RING,SC.ORD.Magnet),1E6*rollMagSupport,'D-','Color',tmpCol(1,:));
	% Fake plot for legend
	pVec(end+1)=plot([-2 -1],[0 0],'-D','Color',tmpCol(1,:));legStr{end+1}='Overall support roll';
	% Plot magnet roll
	pVec(end+1)=plot(findspos(SC.RING,SC.ORD.Magnet),1E6*rollMag,'kx','MarkerSize',8);legStr{end+1}='Magnet roll';


	% Check if girders are registered
	if isfield(SC.ORD,'Girder')
		rollGa=[];rollGb=[];
		% Get girder start rolls
		for ord=SC.ORD.Girder(1,:)
			rollGa(end+1)=SC.RING{ord}.GirderRoll;
		end
		% Get girder ending rolls
		for ord=SC.ORD.Girder(2,:)
			rollGb(end+1)=SC.RING{ord}.GirderRoll;
		end

		% Plot girder roll
		for i=1:size(SC.ORD.Girder,2)
			plot(findspos(SC.RING,SC.ORD.Girder(:,i)),1E6*[rollGa(i),rollGb(i)],'Color',tmpCol(5,:),'LineWidth',4)
		end
		% Fake plot for legend entry
		pVec(end+1)=plot([-2 -1],[0 0],'Color',tmpCol(5,:),'LineWidth',4);legStr{end+1}='Girder roll';
	end

	% Legend and axis stuff
	legend(pVec,legStr);
	set(gca,'xlim',[0 findspos(SC.RING,length(SC.RING)+1)],'box','on')
	ylabel(yLabStr);xlabel('$s$ [m]')

	% Link x-axis
	linkaxes(ax,'x');

	% Make nice
	set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
	set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
	set(findall(gcf,'-property','FontSize'),'FontSize',par.fontSize);
	set(gcf,'color','w');
	drawnow

end
