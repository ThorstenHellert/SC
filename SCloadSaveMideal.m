function [RM1,RM2,RMCO,Bref1,Bref2] = SCloadSaveMideal(SC,flag,varargin)
% SCloadSaveMideal
% ================
%
% NAME
% ----
% SCloadSaveMideal - load/save ideal response matrices and reference trajectory
%
% SYNOPSIS
% --------
% `[RM1,RM2,RMCO,Bref1,Bref2] = SCloadSaveMideal(SC, str)`
%
%
% DESCRIPTION
% -----------
% This script loads the precalculated ideal trajectory response matrices named with the string `str`,
% if they exist. Otherwise the 1 and 2-turn trajectory response and orbit response matrices are
% calculated and saved. Also calculated are the 1 and 2-turn reference trajectories.
%
% INPUTS
% ------
%
% `SC`::
%	SC base structure
% `flag`::
%	String which is used to identify the response matrices
%
% OPTIONS
% -------
% The following options can be specified as name-value pairs:
%
% `'overwrite'` (0)::
%	If true, response matrizes will be calculated and saved even if the files exist
% `'CMords'` (`SC.ORD.CM`):: 
%   List of CM ordinates
% `'BPMords'` (`SC.ORD.BPM`):: 
%   List of BPM ordinates
%
% RETURN VALUE
% ------------
% `RM1`::
% 	1-turn trajectory response matrix.
% `RM2`::
% 	2-turn trajectory response matrix.
% `RMCO`::
% 	Orbit response matrix.
% `Bref1`::
% 	1-turn reference trajectory.
% `Bref2`::
% 	2-turn reference trajectory.
%
% EXAMPLES
% --------
%
% Check if the files `RM1_FODO`, `RM2_FODO` and `RMCO_FODO` exist and load them. Otherwise
% the response matrices for one and two turns and the orbit response matrix are calculated and 
% saved under the previously mentioned file names. Additionally returns the 1 and 2 turn reference 
% trajectory which might be usefull for off-axis injection schemes. 
% ------------------------------------------------------------------
% [RM1,RM2,RMCO,Bref1,Bref2] = SCloadSaveMideal(SC,'_FODO')
% ------------------------------------------------------------------
%
% SEE ALSO
% --------
% *SCgetModelRM*, *SCgetBPMreading*

	
% Parse optional arguments
p = inputParser;
addOptional(p,'overwrite',0);
addOptional(p,'CMords',SC.ORD.CM);
addOptional(p,'BPMords',SC.ORD.BPM);
parse(p,varargin{:});
par=p.Results;

% Do for 1 and 2 turns
for nTurns = 1:2
	% Define file names
	RMfileName = sprintf('RM%d_%s.mat',nTurns,flag);
	RMname     = sprintf('RM%d',nTurns);
	% Check if file can be found
	if exist(RMfileName) == 2 && par.overwrite ~= 1
		% Load response matrix
		load(RMfileName);
		fprintf('Loaded %s from file %s\n',RMname,RMfileName);
	else
		% Calculate response matrix
		eval([genvarname(sprintf('RM%d',nTurns)) '= SCgetModelRM(SC,par.BPMords,par.CMords,''nTurns'',nTurns,''Z0'',SC.INJ.Z0ideal,''useIdealRing'',1);'])
		% Save response matrix into file
		save(RMfileName,genvarname(sprintf('RM%d',nTurns)));
		fprintf('Saved %s into file %s\n',RMname,RMfileName);
	end
end

RMCO = [];
% % % Define file names
% % RMfileName = sprintf('RMCO_%s.mat',flag);
% % RMname     = sprintf('RMCO');
% % if exist(RMfileName) == 2 && par.overwrite ~= 1
% % 	load(RMfileName);
% % 	fprintf('Loaded %s from file %s\n',RMname,RMfileName);
% % else
% % 	% Turn on cavity and radiation
% % 	SC.RING = SCcronoff(SC.RING,'cavityon','radiationon');
% % 	% Calculate response matrix
% % 	eval([genvarname(sprintf('RMCO')) '= SCgetModelRM(SC,par.BPMords,par.CMords,''useIdealRing'',1,''trackMode'',''ORB'');'])
% % 
% % 	save(RMfileName,genvarname(sprintf('RMCO')));
% % 	fprintf('Saved %s into file %s\n',RMname,RMfileName);
% % end

% Calculate 1 and 2 turn reference trajectories
SC.INJ.nParticles = 1;
SC.INJ.nShots     = 1;
SC.INJ.trackMode  = 'TBT';

% Get 1-turn reference trajectory
SC.INJ.nTurns=1;
Bref1 = SCgetBPMreading(SC,'BPMords',par.BPMords);

% Get 2-turn reference trajectory
SC.INJ.nTurns=2;
Bref2 = SCgetBPMreading(SC,'BPMords',par.BPMords);
