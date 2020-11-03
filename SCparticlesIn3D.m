function M = SCparticlesIn3D(R, NPART)
% SCparticlesIn3D
% ===============
%
% NAME
% ----
% SCparticlesIn3D - Re-order particle trajectories
%
%
% DESCRIPTION
% -----------
% `atpass` gives its result in form of an `6x(NPART NELEM*NTURN)` array, where
% `NPART`/`NELEM`/`NTURN` are the number of particles/elements/turns.
% This function takes this output and reorders it into an 3D array with
% dimensions `6x(NELEM*NTURN)xNPART`.
%
% INPUT
% -----
% `R`::
%	Output of `atpass`.
% `NPART`::
%	Number of particles in `R`.
%
% RETURN
% ------
% `M`::
%	3D array in with dimensions `6x(NELEM*NTURN)xNPART`.


nZ = size(R,1);
	NELEM = size(R,2)/NPART;
	M = permute(reshape(R,nZ,NPART,NELEM),[1 3 2]);
end
