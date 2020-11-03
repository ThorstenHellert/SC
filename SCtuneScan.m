function [qSP,SC,maxTurns,finTrans,ERROR] = SCtuneScan(SC,qOrds,qSPvec,varargin)
% SCtuneScan
% ==========
%
% NAME
% ----
% SCtuneScan - Varies quadrupole families to improve beam transmission.
%
% SYNOPSIS
% --------
% `[qSP, SC, maxTurns, finTrans, ERROR] = SCtuneScan(SC, qOrds, qSPvec [, options])`
%
%
% DESCRIPTION
% -----------
% Varies two quadrupole groups specified in cell array `qOrds` on a grid of relative setpoints specified in
% `qSPvec` in a spiral-like pattern to increase the beam transmission. Returns the relative setpoints
% which satisfied the target condition or, if the target could not be reached the values which
% resulted in best transmission.
%
%
% INPUTS
% ------
% `SC`::     `SC` base structure
% `qOrds`::  `[1x2]` cell array of quadrupole ordinates {`[1 x NQ1],[1 x NQ2]`}
% `qSPvec`:: `[1x2]` cell array of quadrupole setpoints {`[SP1_1,...,SP1_N1],[SP2_1,...,SP2_N2]`} with `N2=N1`.
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'nParticles'` (`SC.INJ.nParticles`)::
%   Number of particles used for tracking.
% `'nTurns'` (`SC.INJ.nTurns`)::
%   Number of turns used for tracking.
% `'target'` (`1`)::
%   Transmission target at `'nTurns'`.
% `'fullScan'` (0)::
%	If false, the scan finishes as soon as the target is reached.
% `'plotFlag'` (0)::
%	If true, beam transmission is plotted at every step.
% `'verbose'` (0)::
%	If true, additional information is printed.
%
%
% RETURN VALUES
% -------------
% `qSP`::
%	Final setpoints of quadrupole families (relative to current values)
% `SC`::
%   SC structure with applied setpoints
% `maxTurns`::
%   Array of achieved turns matching the scanning pattern
% `finTrans`::
%   Array of turn-by-turn beam transmission matching the scanning pattern
% `ERROR`::
%	Error value.
%
%
% ERRORS
% ------
% `0`::
% 	Beam transmission target reached.
% `1`::
%	Beam transmission or number of turns increased, target not reached.
% `2`::
% 	Unable to increase beam transmission.
%
%
% SEE ALSO
% --------
% *getBPMReading()*, *SCgenBunches()*


	% Parse optional arguments
	p = inputParser;
	addOptional(p,'verbose',0);
	addOptional(p,'plotFlag',0);
	addOptional(p,'nParticles',SC.INJ.nParticles);
	addOptional(p,'nTurns',SC.INJ.nTurns);
	addOptional(p,'target',1);
	addOptional(p,'fullScan',0);
	parse(p,varargin{:});
	par = p.Results;

	% Input check
	inputCheck();

	% Prelocate output
	maxTurns = nan(length(qSPvec{1}),length(qSPvec{2}));
	finTrans = nan(length(qSPvec{1}),length(qSPvec{2}),par.nTurns);
	ERROR    = 2;
	qSP      = [];
	allInd   = [];


	% Generate the indexing order for the spiral-like search pattern
	tmp           = spiral(max(size(maxTurns)));
	[~,idx]       = sort(tmp(:));
	[q1Ind,q2Ind] = ind2sub([max(size(maxTurns)),max(size(maxTurns))],idx);

	% Start main loop
	for i=1:length(q1Ind)
		q1 = q1Ind(i);
		q2 = q2Ind(i);

		% Set quadrupoles to setpoints
		ords      = horzcat(qOrds{:});
		setpoints = [repmat(qSPvec{1}(q1),1,length(qOrds{1})),repmat(qSPvec{2}(q2),1,length(qOrds{2}))];
		SC = SCsetMags2SetPoints(SC,ords,2,2,setpoints,'method','rel');

		% Calculate beam transmission
		[maxTurns(q1,q2),lostCount,~]= SCgetBeamTransmission(SC,'nParticles',par.nParticles,'nTurns',par.nTurns,'verbose',par.verbose);

		% Store final transmission
		finTrans(q1,q2,:)   = 1-lostCount;

		% Store corresponding setpoint indices
		allInd(end+1,:) = [q1,q2];

		% Plot results
		if p.Results.plotFlag
			plotFunction()
		end

		% If full scan should not be performed
		if ~par.fullScan
			% Check if target is reached and return
			if finTrans(q1,q2,end)>=par.target
				ERROR = 0;
				qSP(1) = qSPvec{1}(q1);
				qSP(2) = qSPvec{2}(q2);
				if par.verbose
					fprintf('Transmission target reached with:\n  %s SetPoint: %.4f\n  %s SetPoint: %.4f\n',SC.RING{qOrds{1}(1)}.FamName,qSP(1),SC.RING{qOrds{2}(1)}.FamName,qSP(2))
				end
				return
			end
		end
	end


	% Check if transmission was improved
	for i=1:size(allInd,1)
		testTrans(i) = finTrans(allInd(i,1),allInd(i,2),end);
		testTurns(i) = maxTurns(allInd(i,1),allInd(i,2));
	end
	% Test final transmission
	[a,b] = sort(testTrans,'descend');
	if a(1)==0
		% Test best number of turns
		[a,b] = sort(testTurns,'descend');
		if a(1)==0
			ERROR=2;
			fprintf('Fail, no transmission at all.\n')
			return
		else
			if par.verbose
				fprintf('No transmission at final turn at all. Best number of turns (%d) reached with:\n  %s SetPoint: %.4f\n  %s SetPoint: %.4f\n',a(1),SC.RING{qOrds{1}(1)}.FamName,qSPvec{1}(allInd(b(1),1)),SC.RING{qOrds{2}(1)}.FamName,qSPvec{2}(allInd(b(1),2)))
			end
		end
	else
		if par.verbose
			fprintf('Transmission target not reached. Best value (%d) reached with:\n  %s SetPoint: %.4f\n  %s SetPoint: %.4f\n',a(1),SC.RING{qOrds{1}(1)}.FamName,qSPvec{1}(allInd(b(1),1)),SC.RING{qOrds{2}(1)}.FamName,qSPvec{2}(allInd(b(1),2)))
		end
	end

	% Get best setpoints
	qSP(1) = qSPvec{1}(allInd(b(1),1));
	qSP(2) = qSPvec{2}(allInd(b(1),2));

	if qSP(1)==qSPvec{1}(q1Ind(1)) && qSP(2)==qSPvec{2}(q2Ind(1))
		fprintf('No improvement possible.\n')
		ERROR = 2;
		return
	else
		ERROR = 1;
	end

	% Set quadrupoles to best setpoints
	ords      = horzcat(qOrds{:});
	setpoints = [repmat(qSP(1),1,length(qOrds{1})),repmat(qSP(2),1,length(qOrds{2}))];
	SC = SCsetMags2SetPoints(SC,ords,2,2,setpoints,'method','rel');

	% Plot scan results
	function plotFunction()
		figure(185);clf
		subplot(2,2,1)
		imagesc(100*finTrans(:,:,end))
		c1=colorbar('northoutside');ylabel(c1,'Beam transmission [%]','Interpreter','none')
		set(gca,'clim',[0 100])%,'ColorScale','log'
		ylabel([SC.RING{qOrds{2}(1)}.FamName ' [rel. to nom setpoint]'],'Interpreter','none');xlabel([SC.RING{qOrds{1}(1)}.FamName ' [rel. to nom setpoint]'],'Interpreter','none')
		tickInd = unique(round([linspace(1,length(qSPvec{1}),5) length(qSPvec{1})]));
		set(gca,'YTick',tickInd,'XTickLabel',qSPvec{1}(tickInd),'YTickLabel',qSPvec{2}(tickInd),'XTick',tickInd)

		subplot(2,2,2)
		imagesc(maxTurns)
		c1=colorbar('northoutside');ylabel(c1,'Number of achieved turns','Interpreter','none')
		set(gca,'ColorScale','log','clim',[1 par.nTurns])
		ylabel([SC.RING{qOrds{2}(1)}.FamName ' [rel. to nom setpoint]'],'Interpreter','none');xlabel([SC.RING{qOrds{1}(1)}.FamName ' [rel. to nom setpoint]'],'Interpreter','none')
		tickInd = unique(round([linspace(1,length(qSPvec{1}),5) length(qSPvec{1})]));
		set(gca,'YTick',tickInd,'XTickLabel',qSPvec{1}(tickInd),'YTickLabel',qSPvec{2}(tickInd),'XTick',tickInd)

		subplot(2,2,[3 4])
		stairs(lostCount);hold on;plot([0 par.nTurns],[SC.INJ.beamLostAt SC.INJ.beamLostAt],'k:')
		set(gca,'xlim',[0 par.nTurns],'ylim',[0 1])
		xlabel('Number of turns');ylabel('EDF of lost count');

		set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
		set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
		set(findall(gcf,'-property','FontSize'),'FontSize',18);set(gcf,'color','w');drawnow
	end

	% Check if input looks reasonable
	function inputCheck()
		if length(qSPvec{1})~=length(qSPvec{2})
			error('Both quad setpoint vectors must be have the same length.\n')
		end
		if par.nTurns==1
			error('Doesn''t work with 1 turn.')
		end
	end
end
