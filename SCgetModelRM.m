function [RM, RING] = SCgetModelRM(SC,BPMords,CMords,varargin)
% SCgetModelRM
% ============
%
% NAME
% ----
% SCgetModelRM - determine the lattice response matrix based on current setpoints
%
% SYNOPSIS
% --------
% `[RM, RING] = SCgetModelRM(SC,BPMords,CMords [, options])`
%
%
% DESCRIPTION
% -----------
% SCgetModelRM calculates the reponse matrix `RM` with the BPMs at the ordinates `BPMords`
% and corrector magnets at ordinates `CMords` using the current magnet setpoints without any
% roll/alignment/calibration errors. `CMords` is a 2-cell array containing the two lists for the
% horizontal and vertical CMs respectively.
% This routine can determine the turn-by-turn RM, as well as the orbit-RM; see
% option 'trackMode'.
%
%
% RETURN VALUES
% -------------
% `RM`::
%	The response matrix given in [m/rad].
% `RING`::
%	The idealised RING structure, which was used to determine the RM.
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'trackMode'` (`'TBT'`)::
%	If `TBT` the turn-by-turn RM is calculated. If `ORB` the orbit RM is calculated, using `findorbit6`
% `'Z0'` (`zeros(6,1)`)::
%	Initial condition for tracking. In `ORB`-mode this is used as the initial guess for `findorbit6`.
% `'nTurns'` (`1`)::
%	Number of turns over which to determine the TBT-RM. Ignored if in `ORB`-mode.
% `'dkick'` (`1e-5`)::
%	Kick [rad] to be added when numerically determining the partial derivatives.
% `'useIdealRing'` (`0`)::
%	If true, the design lattice specified in `SC.IDEALRING` is used.


	% Parse optional arguments
	p = inputParser;
	addOptional(p,'trackMode','TBT');
	addOptional(p,'Z0',zeros(6,1));
	addOptional(p,'nTurns',1);
	addOptional(p,'dkick',1e-5);
	addOptional(p,'useIdealRing',0);
	parse(p,varargin{:});
	trackMode = p.Results.trackMode;
	Z0 = p.Results.Z0(:);
	nTurns = p.Results.nTurns;
	dkick = p.Results.dkick;


	% Check if ideal ring should be used, otherwise get model from current setpoints.
	if p.Results.useIdealRing
		RING = SC.IDEALRING;
	else
		RING = SCgetModelRING(SC);
	end

	% Define used passmethod
	switch trackMode
		case 'TBT'
			trackmethod = @atpass;
			fprintf('Calculating model trajectory response matrix for %d turns.\n',nTurns)
		case 'ORB'
			trackmethod = @orbpass;
			nTurns = 1;
			fprintf('Calculating model orbit response matrix.\n')
		otherwise
			error('trackMode "%s" unknown. Valid values are "TBT" and "ORB".',trackMode)
	end

	% Number of BPMs and CMs
	nBPM = length(BPMords);
	nCM = length([CMords{:}]);

	% Allocate result matrix
	RM = nan(2 * nBPM * nTurns, nCM);
	
	% Calculate initial trajectories
	Ta = trackmethod(RING, Z0, 1, nTurns, BPMords);

	if any(isnan(Ta(:)))
		fprintf('Initial trajectory/orbit is NaN. Aborting. \n')
		return
	end

	PolynomDim={'PolynomB','PolynomA'};
	cnt=1;
	for nDim = 1:2
		for CMord = CMords{nDim}

			if strcmp(RING{CMord}.PassMethod,'CorrectorPass')
				% Backup nominal kick
				KickNominal = RING{CMord}.('KickAngle')(nDim);

				% Vary kick angle
				RING{CMord}.('KickAngle')(nDim) = KickNominal + dkick;

				% Calculate resulting trajectory
				TdB = trackmethod(RING, Z0, 1, nTurns, BPMords);

				% Reset polynom
				RING{CMord}.('KickAngle')(nDim) =  KickNominal;

			else
				% Backup nominal polynom
				PolynomNominal = RING{CMord}.(PolynomDim{nDim});

				% Convert dkick to Polynom-delta
				delta = dkick / RING{CMord}.Length;

				% Vary polynom (positive setpoint -> positive kick -> negative horizontal field)
				RING{CMord}.(PolynomDim{nDim})(1) = PolynomNominal(1)  + (-1)^(nDim) * delta;

				% Calculate resulting trajectory
				TdB = trackmethod(RING, Z0, 1, nTurns, BPMords);

				% Reset polynom
				RING{CMord}.(PolynomDim{nDim}) =  PolynomNominal;

			end

			% Calculate ``derivative''
			dTdB = ( TdB - Ta ) / dkick;

			% Add column to output matrix
			RM(:,cnt) = [dTdB(1,:)'; dTdB(3,:)'];

			cnt=cnt+1;
		end
	end

end

% Auxiliary wrapper to call findorbit6 and atpass with same input arguments
function OUT = orbpass(RING, Z0, newlat, nTurns, REFPTS)
	OUT = findorbit6(RING,REFPTS,Z0);
end
