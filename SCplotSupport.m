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
% *SCregisterSupport*, *SCgetSupportOffset*, *SCgetSupportRoll*, *SCupdateSupport*
	
	% Parse input
	p = inputParser;
	addOptional(p,'fontSize',12);
	addOptional(p,'shiftAxes',0.03);
	addOptional(p,'xLim',[0 findspos(SC.RING,length(SC.RING)+1)]);
	parse(p,varargin{:});
	par = p.Results;
	
	% Check if magnets and BPMs are registered
	if ~isfield(SC.ORD,'Magnet')
		error('Magnets must be registered. Use ''SCregisterMagnets''.')
	elseif ~isfield(SC.ORD,'BPM')
		error('BPMs must be registered. Use ''SCregisterBPMs''.')
	end
	
	% Get s position and support structure offset and rolls with high resolution
	C = findspos(SC.RING,length(SC.RING)+1);
	s = linspace(0,C,100*length(SC.RING));
	offSupportLine  = SCgetSupportOffset(SC,s);
	rollSupportLine = SCgetSupportRoll(SC,s);
	
	% Magnet offsets and rolls
	i=0;
	for ord=SC.ORD.Magnet
		i=i+1;
		% Support structure offset
		offMagSupport(:,i)=SC.RING{ord}.SupportOffset;
		% Support structure roll
		rollMagSupport(:,i)=SC.RING{ord}.SupportRoll;
		
		% Get individual magnet offset
		offMagInd(:,i)=SC.RING{ord}.MagnetOffset;
		% Get individual magnet roll
		rollMagInd(:,i)=SC.RING{ord}.MagnetRoll;
		
		% Get overall magnet offset
		offMagTot(:,i)=SC.RING{ord}.T2([1 3 6]);
		% Get overall magnet roll
		rollMagTot(:,i)=SC.RING{ord}.MagnetRoll + SC.RING{ord}.SupportRoll;
	end
	
	% Loop over individual support structure types
	for type = {'Section','Plinth','Girder'}
		% Check if support structure is registered
		if isfield(SC.ORD,type{1})
			i=1;
			for ordPair=SC.ORD.(type{1})
				% Get girder start and ending offsets
				off.(type{1}).a(:,i)=SC.RING{ordPair(1)}.([type{1} 'Offset']);
				off.(type{1}).b(:,i)=SC.RING{ordPair(2)}.([type{1} 'Offset']);
				% Get girder rolls
				roll.(type{1})(:,i)=SC.RING{ordPair(1)}.([type{1} 'Roll']);
				i = i+1;
			end
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
	yLabOffStr  = {'$\Delta x$ [$\mu$m]','$\Delta y$ [$\mu$m]','$\Delta z$ [$\mu$m]'};
	yLabRollStr = {'$a_z$ [$\mu$rad]','$a_x$ [$\mu$rad]','$a_y$ [$\mu$rad]'};
	titlteOffStr = {'Horizontal Offsets','Vertical Offsets','Longitudinal Offsets'};
	titlteRollStr  = {'Roll (roll around z-axis)','Pitch (roll around x-axis)','Yaw (roll around y-axis)'};
		
	lineSpec.Plinth  = {'Color','r','LineWidth',4};
	lineSpec.Section = {'Color',tmpCol(7,:),'LineWidth',2,'LineStyle',':'};
	lineSpec.Girder  = {'Color',tmpCol(5,:),'LineWidth',4};

	% Loop over both transverse planes
	for nDim=1:3
		

		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot individual offset contributions
		ax(3*(nDim-1)+1,1)=subplot(12,2,2*4*(nDim-1)+ [1 3]);hold on;
		pVec=[];legStr=[];
		
		% Plot line of support offset
		stairs(s,1E6*offSupportLine(nDim,:));
		% Plot support offset at magnets
		plot(findspos(SC.RING,SC.ORD.Magnet),1E6*offMagSupport(nDim,:),'D','Color',tmpCol(1,:));
		% Fake plot for legend
		pVec(end+1)=plot([-2 -1],[0 0],'-D','Color',tmpCol(1,:));legStr{end+1}='Overall support structure';
		
		% Plot individual magnet offset
		pVec(end+1)=plot(findspos(SC.RING,SC.ORD.Magnet),1E6*offMagInd(nDim,:),'kx','MarkerSize',8);legStr{end+1}='Individual Magnet';
		
		% Loop over support structure types
		for type = {'Section','Plinth','Girder'}
			% Check if support structure is registered
			if isfield(SC.ORD,type{1})
				% Plot plinth offset
				for i=1:size(SC.ORD.(type{1}),2)
					% Check if support structure spans over injection point
					if diff(findspos(SC.RING,SC.ORD.(type{1})(:,i)))<0
						for nCase=1:2
							if nCase==1
								% Interpolate between last support structure and end of ring
								splot  = findspos(SC.RING,[SC.ORD.(type{1})(1,i) length(SC.RING)]);
								sint   = [findspos(SC.RING,SC.ORD.(type{1})(1,i)),findspos(SC.RING,SC.ORD.(type{1})(2,i))+C];
							else
								% Interpolate between beginning of ring and 1st support structure
								splot  = findspos(SC.RING,[1 SC.ORD.(type{1})(2,i)]);
								sint   = [-findspos(SC.RING,SC.ORD.(type{1})(2,i)),findspos(SC.RING,SC.ORD.(type{1})(2,i))];
							end
							offInt = interp1(sint,[off.(type{1}).a(nDim,i) off.(type{1}).b(nDim,i)],splot);
							plot(splot,1E6*offInt,lineSpec.(type{1}){:})
						end
					else
						plot(findspos(SC.RING,SC.ORD.(type{1})(:,i)),1E6*[off.(type{1}).a(nDim,i) off.(type{1}).b(nDim,i)],lineSpec.(type{1}){:})
					end
				end
				% Fake plot for legend entry
				pVec(end+1)=plot([-2 -1],[0 0],lineSpec.(type{1}){:});legStr{end+1}=sprintf('Individual %s',type{1});
			end
		end
		
		% Legend and axis stuff
		legend(pVec,legStr);
		set(gca,'xlim',par.xLim,'box','on')%,'XTickLabel','')
		ylabel(yLabOffStr{nDim});
		title(titlteOffStr{nDim})
		
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot overall magnet offset
		ax(3*(nDim-1)+2,1)=subplot(12,2,2*4*(nDim-1)+ 5);
		
		% Plot magnet offset
		plot(findspos(SC.RING,SC.ORD.Magnet),1E6*offMagTot(nDim,:),'kO-');
		
		% Legend and axis stuff
		legend('Overall magnet offset');
		set(gca,'xlim',par.xLim,'box','on')%,'XTickLabel','')
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
		set(gca,'xlim',par.xLim,'box','on')%,'XTickLabel','')
		ylabel(yLabOffStr{nDim});
		if nDim==3
			xlabel('$s$ [m]')
		else
% 			set(gca,'XTickLabel','')
		end

		
		
		
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot individual roll contributions
		ax(3*(nDim-1)+1,2)=subplot(12,2,2*4*(nDim-1)+ [2 4]);hold on;
		pVec=[];legStr=[];
					
		% Plot line of total support rolls
		stairs(s,1E6*rollSupportLine(nDim,:),'Color',tmpCol(1,:));
		% Plot total support rolls at magnets
		plot(findspos(SC.RING,SC.ORD.Magnet),1E6*rollMagSupport(nDim,:),'D','Color',tmpCol(1,:));
		% Fake plot for legend
		pVec(end+1)=plot([-2 -1],[0 0],'-D','Color',tmpCol(1,:));legStr{end+1}='Overall support structure';
		
		% Plot magnet individual roll
		pVec(end+1)=plot(findspos(SC.RING,SC.ORD.Magnet),1E6*rollMagInd(nDim,:),'kx','MarkerSize',8);legStr{end+1}='Individual Magnet';
		
		% Loop over support structure types
		for type = {'Section','Plinth','Girder'}
			% Check if support structure is registered
			if isfield(SC.ORD,type{1})
				% Plot plinth offset
				for i=1:size(SC.ORD.(type{1}),2)
					% Check if support structure spans over injection point
					if diff(findspos(SC.RING,SC.ORD.(type{1})(:,i)))<0
						for nCase=1:2
							if nCase==1
								% Interpolate between last support structure and end of ring
								splot  = findspos(SC.RING,[SC.ORD.(type{1})(1,i) length(SC.RING)]);
								sint   = [findspos(SC.RING,SC.ORD.(type{1})(1,i)),findspos(SC.RING,SC.ORD.(type{1})(2,i))+C];
							else
								% Interpolate between beginning of ring and 1st support structure
								splot  = findspos(SC.RING,[1 SC.ORD.(type{1})(2,i)]);
								sint   = [-findspos(SC.RING,SC.ORD.(type{1})(2,i)),findspos(SC.RING,SC.ORD.(type{1})(2,i))];
							end
							rollInt = interp1(sint,roll.(type{1})(nDim,i)*[1 1],splot);
							plot(splot,1E6*rollInt,lineSpec.(type{1}){:})
						end
					else
						plot(findspos(SC.RING,SC.ORD.(type{1})(:,i)),1E6*roll.(type{1})(nDim,i)*[1 1],lineSpec.(type{1}){:})
					end
				end
				% Fake plot for legend entry
				pVec(end+1)=plot([-2 -1],[0 0],lineSpec.(type{1}){:});legStr{end+1}=sprintf('Individual %s',type{1});
			end
		end
		

		
		% Legend and axis stuff
		legend(pVec,legStr);
		set(gca,'xlim',par.xLim,'box','on','YAxisLocation','right')%,'XTickLabel',''
		ylabel(yLabRollStr{nDim});
		title(titlteRollStr{nDim})
		
				
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plot overall magnet roll
		ax(3*(nDim-1)+2,2)=subplot(12,2,2*4*(nDim-1)+ 6);hold on
				
		% Plot magnet roll
		plot(findspos(SC.RING,SC.ORD.Magnet),1E6*rollMagTot(nDim,:),'kO-');
		
		% Legend and axis stuff
		legend('Overall magnet roll');
		set(gca,'xlim',par.xLim,'box','on','YAxisLocation','right')%,'XTickLabel',''
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
		set(gca,'xlim',par.xLim,'box','on','YAxisLocation','right')
		ylabel(yLabRollStr{nDim});
		if nDim==3
			xlabel('$s$ [m]')
		end
		
	end
	
	% Link x-axis
	linkaxes(ax,'x');
	
	% Rearrange plots (not sure if it works on your screen resolution)
	for nDim=1:3
		for nAx=1:3
			for n=1:2
				set(ax(nAx+3*(nDim-1),n),'Position',get(ax(nAx+3*(nDim-1),n),'Position') - ((nDim-1)-0.4*(nAx-1))*[0 par.shiftAxes 0 0])
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
