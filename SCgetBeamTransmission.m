function [maxTurns,lostCount,ERROR] = SCgetBeamTransmission(SC,varargin)
% SCgetBeamTransmission
% =====================
%
% NAME
% ----
% SCgetBeamTransmission - Calculate turn-by-turn beam transmission
%
% SYNOPSIS
% --------
% `[maxTurns, lostCount, ERROR] = SCgetBeamTransmission(SC [, options])`
%
%
% DESCRIPTION
% -----------
% Calculates the turn-by-turn beam transmission. A bunch is generated according
% to the injection setting and tracking is performed. At each turn the relative
% amount of lost particles is calculated. The number of survived turns is
% determined by using the user specified field `SC.INJ.beamLostAt`.
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'nParticles'` (`SC.INJ.nParticles`)::
%	Number of particles.
% `'nTurns'` (`SC.INJ.nTurns`)::
% 	Number of turns.
% `'plotFlag'` (0)::
%	If true, particle lost count is plotted.
% `'verbose'` (0)::
%	If true, additional information is printed.
%
%
% RETURN VALUE
% ------------
% `maxTurns`::
% 	Number of achieved turns.
% `lostCount`::
% 	Turn-by-turn particle losses.
% `ERROR`::
% 	False if beam survives all turns.
%
% SEE ALSO
% --------
% *SCgenBunches*, *SCgetBPMreading*


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse optional arguments
p = inputParser;
addOptional(p,'nParticles',SC.INJ.nParticles);
addOptional(p,'nTurns',SC.INJ.nTurns);
addOptional(p,'plotFlag',0);
addOptional(p,'verbose',0);
parse(p,varargin{:});
par = p.Results;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization 

ERROR = 1;

% Change injection setting
SC.INJ.nParticles = par.nParticles;
SC.INJ.nTurns	  = par.nTurns;

% Print information
if par.verbose
	fprintf('Calculating maximum beam transmission for %d particles and %d turns: ',SC.INJ.nParticles,SC.INJ.nTurns)
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function 

% Generate bunches
Zin = SCgenBunches(SC);

% Tracking
T = atpass(SC.RING, Zin,1,SC.INJ.nTurns,length(SC.RING)+1);

% Check for multiparticle tracking
if SC.INJ.nParticles > 1

	% Reshape trajectories in tensor
	M = SCparticlesIn3D(T,SC.INJ.nParticles);

	% Read x position
	Tx = squeeze(M(1,:,:));

	% Check for detactable BPM signal
	maxTurns = find( sum( isnan( Tx ),2 ) > ( SC.INJ.nParticles * SC.INJ.beamLostAt ),1)-1;

else
	% Read x position
	Tx = T(1,:)';

	% Check for detactable BPM signal
	maxTurns = find( isnan( Tx ) ,1)-1;
end

% Set to max number if beam is not lost
if isempty(maxTurns)
	maxTurns = SC.INJ.nTurns;
	ERROR    = 0;
end

% Calculate statistics
lostCount = sum( isnan( Tx ),2 )/SC.INJ.nParticles;

% Plot particle lost count
if par.plotFlag
	figure(12),clf
	stairs(lostCount);hold on;plot([0 SC.INJ.nTurns],[SC.INJ.beamLostAt SC.INJ.beamLostAt],'k:')
	set(gca,'xlim',[0 SC.INJ.nTurns],'ylim',[0 1])
	xlabel('Number of turns');ylabel('EDF of lost count');
	set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
	set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
	set(findall(gcf,'-property','FontSize'),'FontSize',18);
	set(gcf,'color','w');
	drawnow
end
% Print information
if par.verbose
	fprintf('%d turns and %.0f%% transmission.\n',maxTurns,100*(1-lostCount(end)))
end


end


