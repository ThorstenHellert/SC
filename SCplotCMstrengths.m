function  SCplotCMstrengths(SC,varargin)
% SCplotCMstrengths
% =================
%
% NAME
% ----
% SCplotCMstrengths - plots the current CM setpoints
%
% SYNOPSIS
% --------
% `SCplotCMstrengths(SC)`
%
%
% DESCRIPTION
% -----------
% Plots the current CM setpoints.
%
% INPUTS
% ------
% `SC`:: SC base structure
%
% SEE ALSO
% --------
% *SCsetCMs2SetPoints*


	% Define Hor/Ver CM setpoint field name
	fieldNames = {'SetPointB','SetPointA'}; % {hor,ver}

	% Create figure
	figure(86);clf

	% Loop over transverse planes
	for nDim=1:2

		% Get CM s-position
		CMs{nDim} = findspos(SC.RING,SC.ORD.CM{nDim})';

		% Get CM strength (kick angle or PolynomA/B times lengthof magnet)
		nCM=1;
		for ord=SC.ORD.CM{nDim}
			if strcmp(SC.RING{ord}.PassMethod,'CorrectorPass')
				CMval{nDim}(nCM) = SC.RING{ord}.KickAngle(nDim);
			else
				CMval{nDim}(nCM) = SC.RING{ord}.(fieldNames{nDim})(1) .* SC.RING{ord}.Length;
			end
			nCM=nCM+1;
		end
	end

	% Plot CM strength along the ring % % % % % % % % % % % % % % % % % % % % % % % % % % %
	subplot(2,1,1);hold on
	bar(CMs{1},1E6*CMval{1});
	bar(CMs{2},1E6*CMval{2});

	ylabel('CM strengh [$\mu$rad]');xlabel('s [m]');
	title([sprintf('$\\Theta_x=%.0f\\mu$rad rms,            ', 1E6*sqrt(mean(CMval{1}.^2))),...
		sprintf('$\\Theta_y=%.0f\\mu$rad rms', 1E6*sqrt(mean(CMval{2}.^2)))],'Interpreter','latex')
	set(gca,'box','on')

	% Plot CDF of all CM strengths % % % % % % % % % % % % % % % % % % % % % % % % % % %
	subplot(2,1,2);hold on
	[ah,bh]=hist(1E6*(CMval{1}),length(CMval{1}));
	[av,bv]=hist(1E6*(CMval{2}),length(CMval{2}));
	stairs(bh,cumsum(ah)./sum(ah),'LineWidth',2);
	stairs(bv,cumsum(av)./sum(av),'LineWidth',2);

	% Labels and limits
	xlabel('CM strengh [$\mu$rad]');ylabel('CDF');
	legend({'Horizontal','Vertical'})
	set(gca,'ylim',[0 1],'box','on')

	% Make chic
	set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
	set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
	set(findall(gcf,'-property','FontSize'),'FontSize',18);
	set(gcf,'color','w');
	drawnow

