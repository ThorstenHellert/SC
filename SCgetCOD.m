function COD = SCgetCOD(SC,varargin)
% SCgetCOD
% ========
%
% NAME
% ----
% SCgetCOD - Calculates closed orbit deviations with respect to magnet centers.
%
% SYNOPSIS
% --------
% `COD = SCgetCOD(SC [, options])`
%
% DESCRIPTION
% -----------
% This function calculates the closed orbit deviation from the magnetic centers.
%
% INPUTS
% ------
%
% `SC`::
%	SC base structure
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'ords'` (`SC.ORD.Magnet`):: 
% List of magnet ordinates at which the orbit deviation should be evaluated.
%
% RETURN VALUES
% -------------
% `COD`::
%	Orbit deviation w.r.t. magnet centers
%
% EXAMPLES
% --------
% Calculate closed orbit deviation in magnets named 'SF' and 'SD' and plot results.
% ------------------------------------------------------------------
% COD = SCgetCOD(SC,'ords',SCgetOrds(SC.RING,'SF|SD'),'plot',1);
% ------------------------------------------------------------------

	
	
	% Parse input
	p = inputParser;
	addOptional(p,'ords',SC.ORD.Magnet);
	addOptional(p,'plot',0);
	parse(p,varargin{:});
	par = p.Results;
	
	
	% Calculate orbit
	T = findorbit6(SC.RING,par.ords);
	if any(isnan(T))
		warning('Closed orbit could not be found.')
		COD = nan(2,length(par.ords));
		return
	end

	% Get magnet offsets
	magOffset=zeros(2,length(par.ords));i=1;
	for ord=par.ords
		magOffset(:,i) = SC.RING{ord}.T2([1 3]);
		i=i+1;
	end
	
	% Calculate orbit deviation in magnets
	COD = T([1,3],:) - magOffset;
	

	% Plot results?
	if par.plot	
		
		% Get s-positions
		sPos = findspos(SC.RING,par.ords);
		
		ylabelStr = {'$\Delta x$ [mm]','$\Delta y$ [mm]'};
		figure(784);clf
		for nDim=1:2
			ax(nDim) = subplot(2,1,nDim);hold on
			plot(sPos,1E3*COD(nDim,:),'LineWidth',2)
			plot(sPos,1E3*magOffset(nDim,:),'LineWidth',2)
			plot(sPos,1E3*T(2*nDim-1,:),'LineWidth',2)
			
			legend({'COD in magnets','Magnet offset','Orbit'})
			xlabel('s [m]');
			ylabel(ylabelStr{nDim});
			set(gca,'box','on','xlim',[min(sPos),max(sPos)])
		end
		set(findall(gcf,'-property','FontSize'),'FontSize',18);
		set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
		set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
		set(gcf,'color','w');
		drawnow
		
		% Link x-axis
		linkaxes([ax(1) ax(2)],'x')
	end
	
end

