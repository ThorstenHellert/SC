function Z = SCgenBunches(SC)
% SCgenBunches
% ============
%
% NAME
% ----
% SCgenBunches - Generate injected beam particle coordinates.
%
% SYNOPSIS
% --------
% `Z = SCgenBunches(SC)`
%
%
% DESCRIPTION
% -----------
% Generates bunches according to the current injection setup as defined in
% `SC.INJ`. The injection setup includes the fields:
%
% `nParticles`::
%	Number of particles for tracking.
% `nTurns`::
%	Number of turns for tracking.
% `nShots`::
%	Number of injections used for averaging the BPM readings.
% `Z0ideal`::
%	[6 x 1] array defining the ideal injection vector.
% `Z0`::
%	[6 x 1] array defining the current injection vector.
% `beamSize`::
%	[6 x 6] array defining the beam sigma matrix.
% `randomInjectionZ`::
%	[6 x 1] array defining the shot-to-shot injection jitter.
% `trackMode`::
%	String defining the tracking mode.
%
% For each bunch the random injection error is calculated and added to the mean
% injected beam trajectory. If the number of particles per bunch is larger than
% one, individual particle launch conditions are randomly distributed around
% the bunch centroid using the beam sigma matrix.  Otherwise the single
% macro-particle per bunch is launched at the bunch centroid trajectory. The
% global injection count `SCinjections` is increased appropriately. If a
% post-injection function is defined in `SC.INJ.postFun`, it is applied on the
% generated particles.
%
% INPUTS
% ------
%
% `SC`::
%	SC base structure
%
%
% RETURN VALUES
% -------------
% `Z`::
%	[6 x (#Bunches x #Shots x #Particles)] array of injected beam particle coordinates
%
%
% SEE ALSO
% --------
% *SCgetBPMreading*


	global SCinjections

	% Calculate bunch centroid trajectory
	Z = repmat(SC.INJ.randomInjectionZ.*SCrandnc(2,6,1) + SC.INJ.Z0 ,1,SC.INJ.nParticles);

    % If not only single macro particle is used
    if SC.INJ.nParticles~=1
		% Perform eigenvalue decomposition -> chol() doesn't work if sigma is not positive definite
		[V,L] = eig(SC.INJ.beamSize);
		% Generate particles from sigma matrix
		particles = V*sqrt(L) * SCrandnc(3,6,SC.INJ.nParticles);
		% Add to bunch centroid coordinates
		Z = Z + particles;
	end

	% Apply post injection function, if defined
	if isfield(SC.INJ,'postFun') && isa(SC.INJ.postFun, 'function_handle')
		Z = SC.INJ.postFun(Z);
	end

	% Increase total number of injected beams
	SCinjections = SCinjections + 1;

end

