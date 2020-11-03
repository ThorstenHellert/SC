function SC = SCinit(RING)
% SCinit
% ======
%
% NAME
% ----
% SCinit - Initializes the toolbox
%
% SYNOPSIS
% --------
% `SC = SCinit(RING)`
%
% DESCRIPTION
% -----------
% This function is used to initialize the toolbox. Required input `RING` is an AT lattice cell structure,
% which is also saved as the ideal lattice in `SC.IDEALRING`. The global variables `plotFunctionFlag`
% and `SCinjections` are set to initial values.
%
% GLOBALS
% -------
% `plotFunctionFlag` ([])::
% 	If not empty, every call of *SCgetBPMreading* calls *SCplotBPMreading*, hence plotting particle
%   trajectories and BPM readings.
% `SCinjections` (0)::
% 	Counts the number of injected beams. Gets increased by *SCgenBunches*
%
% RETURN VALUE
% ------------
% `SC`::
% 	Base structure including the lattice structure.


	global plotFunctionFlag SCinjections

	% Store lattice
	SC.RING = RING;

	% Define ideal lattice
	SC.IDEALRING = RING;

	% Relative amount of partcles which may be lost before BPM reading is NaN
	SC.INJ.beamLostAt = 1;

	% Design injected trajectory
	SC.INJ.Z0ideal = zeros(6,1);

	% Injetced trajectory
	SC.INJ.Z0 = SC.INJ.Z0ideal;

	% Injected bunch beam size
	SC.INJ.beamSize = zeros(6,6);

	% Injected beam random trajectory jitter
	SC.INJ.randomInjectionZ = zeros(6,1);

	% Number of particles per bunch
	SC.INJ.nParticles = 1;

	% Number of turns for tracking
	SC.INJ.nTurns = 1;

	% Number of injections for averaging BPM reading
	SC.INJ.nShots = 1;

	% Tracking mode
	SC.INJ.trackMode = 'TBT';

	% Total number of injeceted beams
	SCinjections = 0;

	% Set plot function flag
	plotFunctionFlag = [];	

end
