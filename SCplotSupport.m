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
% This function plots the overall offsets [dx,dy,dz] and rolls [az,ax,ay] of all magnets and BPMs,
% as well as the individual contributions from different support structures (if registered).
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
% `'fontSize'` (12)::
%   Figure font size.
% `'shiftAxes'` (0.03)::
%   Axes are reanranged for grouping. Depending on screen resolution this value may be adjusted.
%
%
% SEE ALSO
% --------
% *SCgetSupportOffset*, *SCupdateSupport*
	
	% Parse input
	p = inputParser;
	addOptional(p,'fontSize',12);
	addOptional(p,'shiftAxes',0.03);
	parse(p,varargin{:});
	par = p.Results;
	
	% Check if magnets and BPMs are registered
	if ~isfield(SC.ORD,'Magnet')
		error('Magnets must be registered. Use ''SCregisterMagnets''.')
	elseif ~isfield(SC.ORD,'BPM')
		error('BPMs must be registered. Use ''SCregisterBPMs''.')
	end
	
	
	
	% Get s position and support structure offset
	C = findspos(SC.RING,length(SC.RING)+1);
	s = linspace(0,C,100*length(SC.RING));
	off = SCgetSupportOffset(SC,s);
	% Magnet offsets and rolls
	i=0;
	for ord=SC.ORD.Magnet
		i=i+1;
		% Support structure offset
		offMagSupport(1:3,i)=SC.RING{ord}.SupportOffset(1:3);
		% Support structure roll
		rollMagSupport(1:3,i)=SC.RING{ord}.SupportRoll(1:3);
		
		% Get magnet-to-magnet offset
		offMagInd(1:3,i)=SC.RING{ord}.MagnetOffset(1:3);
		% Get magnet-to-magnet roll
		rollMagInd(1:3,i)=SC.RING{ord}.MagnetRoll(1:3);
		
		% Get overall magnet offset
		offMag(1:3,i)=SC.RING{ord}.T2([1 3 6]);
		% Get overall magnet roll
		rollMag(1:3,i)=SC.RING{ord}.MagnetRoll + SC.RING{ord}.SupportRoll;

	end
	
	% Check if sections are registered
	if isfield(SC.ORD,'Section')
		i=1;
		% Get section offsets
		for ordPair=SC.ORD.Section
			offSa(1:3,i)=SC.RING{ordPair(1)}.SectionOffset(1:3);
			offSb(1:3,i)=SC.RING{ordPair(2)}.SectionOffset(1:3);
			i = i+1;
		end
	end
	if isfield(SC.ORD,'Plinth')
		i=1;
		% Get plinth  offsets
		for ordPair=SC.ORD.Plinth
			offPa(1:3,i)=SC.RING{ordPair(1)}.PlinthOffset(1:3);
			offPb(1:3,i)=SC.RING{ordPair(2)}.PlinthOffset(1:3);
			i = i+1;
		end
	end
	% Check if girders are registered
	if isfield(SC.ORD,'Girder')
		i=1;
		for ordPair=SC.ORD.Girder
			% Get girder start and ending offsets
			offGa(1:3,i)=SC.RING{ordPair(1)}.GirderOffset(1:3);
			offGb(1:3,i)=SC.RING{ordPair(2)}.GirderOffset(1:3);
			% Get girder rolls
			rollGa(1:3,i)=SC.RING{ordPair(1)}.GirderRoll(1:3);
			i = i+1;
		end
	end
	% BPM offsets
	sBPM         = findspos(SC.RING,SC.ORD.BPM);
	offBPM       = cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'Offset'));
	offBPMStruct = cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'SupportOffset'));
	offBPM(      :,3) = 0;
	offBPMStruct(:,3) = 0;
	% BPM rolls
	rollBPM       = atgetfieldvalues(SC.RING(SC.ORD.BPM),'Roll',{1,1});
	rollBPMStruct = atgetfieldvalues(SC.RING(SC.ORD.BPM),'SupportRoll',{1,1});
	rollBPM(      :,2:3) = 0;
	rollBPMStruct(:,2:3) = 0;
	
	
	% Create figure
	figure(1213);clf;tmpCol=get(gca, 'ColorOrder');ax=[];
	yLabOffStr = {'$\Delta x$ [$\mu$m]','$\Delta y$ [$\mu$m]','$\Delta z$ [$\mu$m]'};
	yLabRollStr = {'$a_z$ [$\mu$rad]','$a_x$ [$\mu$rad]','$a_y$ [$\mu$rad]'};
	titleStr = {'Horizontal','Vertical','Longitudinal'};
	rollStr = {'Roll (roll around z-axis)','Pitch (roll around x-axis)','Yaw (roll around y-axis)'};

	% Loop over both transverse planes
	for nDim=1:3
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Offsets
		
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot individual offset contributions
		ax(3*(nDim-1)+1,1)=subplot(12,2,2*4*(nDim-1)+ [1 3]);hold on;
		pVec=[];legStr=[];
		
		
		
		% Plot line of support offset
		stairs(s,1E6*off(nDim,:));
			
		% Plot support offset
		plot(findspos(SC.RING,SC.ORD.Magnet),1E6*offMagSupport(nDim,:),'D','Color',tmpCol(1,:));
		% Fake plot for legend
		pVec(end+1)=plot([-2 -1],[0 0],'-D','Color',tmpCol(1,:));legStr{end+1}='Overall support structure';
		% Plot magnet offset
		pVec(end+1)=plot(findspos(SC.RING,SC.ORD.Magnet),1E6*offMagInd(nDim,:),'kx','MarkerSize',8);legStr{end+1}='Individual Magnet';
		

		% Check if plinths are registered
		if isfield(SC.ORD,'Plinth')
			% Plot plinth offset
			for i=1:size(SC.ORD.Plinth,2)
				% Check if support structure spans over injection point
				if diff(findspos(SC.RING,SC.ORD.Plinth(:,i)))<0
					for nCase=1:2
						if nCase==1
							% Interpolate between last support structure and end of ring
							splot  = findspos(SC.RING,[SC.ORD.Plinth(1,i) length(SC.RING)]);
							sint   = [findspos(SC.RING,SC.ORD.Plinth(1,i)),findspos(SC.RING,SC.ORD.Plinth(2,i))+C];
						else
							% Interpolate between beginning of ring and 1st support structure
							splot  = findspos(SC.RING,[1 SC.ORD.Plinth(2,i)]);
							sint   = [-findspos(SC.RING,SC.ORD.Plinth(2,i)),findspos(SC.RING,SC.ORD.Plinth(2,i))];
						end
						offInt = interp1(sint,[offPa(nDim,i) offPb(nDim,i)],splot);
						plot(splot,1E6*offInt,'Color','r','LineWidth',4)
					end
				else
					plot(findspos(SC.RING,SC.ORD.Plinth(:,i)),1E6*[offPa(nDim,i) offPb(nDim,i)],'Color','r','LineWidth',4)
				end
			end
			% Fake plot for legend entry
			pVec(end+1)=plot([-2 -1],[0 0],'Color','r','LineWidth',4);legStr{end+1}='Individual Plinth';
		end
		% Check if girders are registered
		if isfield(SC.ORD,'Girder')
			% Plot girder offset
			for i=1:size(SC.ORD.Girder,2)
				plot(findspos(SC.RING,SC.ORD.Girder(:,i)),1E6*[offGa(nDim,i),offGb(nDim,i)],'Color',tmpCol(5,:),'LineWidth',4)
			end
			% Fake plot for legend entry
			pVec(end+1)=plot([-2 -1],[0 0],'Color',tmpCol(5,:),'LineWidth',4);legStr{end+1}='Individual Girder';
		end
		% Check if sections are registered
		if isfield(SC.ORD,'Section')
			% Plot section offset
			for i=1:size(SC.ORD.Section,2)
				% Check if support structure spans over injection point
				if diff(findspos(SC.RING,SC.ORD.Section(:,i)))<0
					for nCase=1:2
						if nCase==1
							% Interpolate between last support structure and end of ring
							splot  = findspos(SC.RING,[SC.ORD.Section(1,i) length(SC.RING)]);
							sint   = [findspos(SC.RING,SC.ORD.Section(1,i)),findspos(SC.RING,SC.ORD.Section(2,i))+C];
						else
							% Interpolate between beginning of ring and 1st support structure
							splot  = findspos(SC.RING,[1 SC.ORD.Section(2,i)]);
							sint   = [-findspos(SC.RING,SC.ORD.Section(2,i)),findspos(SC.RING,SC.ORD.Section(2,i))];
						end
						offInt = interp1(sint,[offSa(nDim,i) offSb(nDim,i)],splot);
						plot(splot,1E6*offInt,'Color',tmpCol(7,:),'LineWidth',2,'LineStyle',':')
					end
				else
					plot(findspos(SC.RING,SC.ORD.Section(:,i)),1E6*[offSa(nDim,i) offSb(nDim,i)],'Color',tmpCol(7,:),'LineWidth',2,'LineStyle',':')
				end
			end
			% Fake plot for legend entry
			pVec(end+1)=plot([-2 -1],[0 0],'Color',tmpCol(7,:),'LineWidth',3,'LineStyle',':');legStr{end+1}='Individual Section';
		end
		
		% Legend and axis stuff
		legend(pVec,legStr);
		set(gca,'xlim',[0 findspos(SC.RING,length(SC.RING)+1)],'box','on','XTickLabel','')
		ylabel(yLabOffStr{nDim});
		title(sprintf('%s Offsets',titleStr{nDim}))
		
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot overall magnet offset
		ax(3*(nDim-1)+2,1)=subplot(12,2,2*4*(nDim-1)+ 5);
		
		% Plot magnet offset
		plot(findspos(SC.RING,SC.ORD.Magnet),1E6*offMag(nDim,:),'kO-');
		
		% Legend and axis stuff
		legend('Overall magnet offset');
		set(gca,'xlim',[0 findspos(SC.RING,length(SC.RING)+1)],'box','on','XTickLabel','')
		ylabel(yLabOffStr{nDim});
		
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot BPM offset
		ax(3*(nDim-1)+3,1)=subplot(12,2,2*4*(nDim-1)+ 7);hold on
		
		% Plot random BPM offset
		plot(sBPM,1E6*offBPM(:,nDim),'O','Color',tmpCol(2,:),'MarkerSize',6);
		% Plot BPM offset from support structure
		plot(sBPM,1E6*offBPMStruct(:,nDim),'-','Color',tmpCol(2,:));
		
		% Legend and axis stuff
		legend({'Random BPM offset','BPM support offset'});
		set(gca,'xlim',[0 findspos(SC.RING,length(SC.RING)+1)],'box','on')%,'XTickLabel','')
		ylabel(yLabOffStr{nDim});
		if nDim==3
			xlabel('$s$ [m]')
		else
			set(gca,'XTickLabel','')
		end

		
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Rolls
		
			
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot individual roll contributions
		ax(3*(nDim-1)+1,2)=subplot(12,2,2*4*(nDim-1)+ [2 4]);hold on;
		pVec=[];legStr=[];
		
				
		% Plot support roll
		plot(findspos(SC.RING,SC.ORD.Magnet),1E6*rollMagSupport(nDim,:),'D-','Color',tmpCol(1,:));
		% Fake plot for legend
		pVec(end+1)=plot([-2 -1],[0 0],'-D','Color',tmpCol(1,:));legStr{end+1}='Overall support structure';
		% Plot magnet roll
		pVec(end+1)=plot(findspos(SC.RING,SC.ORD.Magnet),1E6*rollMagInd(nDim,:),'kx','MarkerSize',8);legStr{end+1}='Individual Magnet';
		
		% Check if girders are registered
		if isfield(SC.ORD,'Girder')
			% Plot girder roll
			for i=1:size(SC.ORD.Girder,2)
				plot(findspos(SC.RING,SC.ORD.Girder(:,i)),1E6*rollGa(nDim,i)*[1 1],'Color',tmpCol(5,:),'LineWidth',4)
			end
			% Fake plot for legend entry
			pVec(end+1)=plot([-2 -1],[0 0],'Color',tmpCol(5,:),'LineWidth',4);legStr{end+1}='Individual Girder';
		end
		
		% Legend and axis stuff
		legend(pVec,legStr);
		set(gca,'xlim',[0 findspos(SC.RING,length(SC.RING)+1)],'box','on','XTickLabel','','YAxisLocation','right')
		ylabel(yLabRollStr{nDim});
		title(rollStr{nDim})
		
				
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot overall magnet roll
		ax(3*(nDim-1)+2,2)=subplot(12,2,2*4*(nDim-1)+ 6);hold on
				
		% Plot magnet roll
		plot(findspos(SC.RING,SC.ORD.Magnet),1E6*rollMag(nDim,:),'kO-');
		
		% Legend and axis stuff
		legend('Overall magnet roll');
		set(gca,'xlim',[0 findspos(SC.RING,length(SC.RING)+1)],'box','on','XTickLabel','','YAxisLocation','right')
		ylabel(yLabRollStr{nDim});
		
			
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot BPM Roll
		ax(3*(nDim-1)+3,2)=subplot(12,2,2*4*(nDim-1)+ 8);hold on
		
		% Plot random BPM roll
		plot(sBPM,1E6*rollBPM(:,nDim),'O','Color',tmpCol(2,:),'MarkerSize',6);
		% Plot BPM roll from support structure
		plot(sBPM,1E6*rollBPMStruct(:,nDim),'-','Color',tmpCol(2,:));
		
		% Legend and axis stuff
		legend({'Random BPM roll','BPM support roll'});
		set(gca,'xlim',[0 findspos(SC.RING,length(SC.RING)+1)],'box','on','YAxisLocation','right')
		ylabel(yLabRollStr{nDim});
		if nDim==3
			xlabel('$s$ [m]')
		else
			set(gca,'XTickLabel','')
		end
	end
	
	% Link x-axis
	linkaxes(ax,'x');
	
	% Rearrange plots
	for nDim=1:3
		for nAx=1:3
			for n=1:2
				set(ax(nAx+3*(nDim-1),n),'Position',get(ax(nAx+3*(nDim-1),n),'Position') - (nDim-1)*[0 par.shiftAxes 0 0])
			end
		end
	end
	
	% Make nice
	set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
	set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
	set(findall(gcf,'-property','FontSize'),'FontSize',par.fontSize);
	set(gcf,'color','w');
	drawnow
	
end
