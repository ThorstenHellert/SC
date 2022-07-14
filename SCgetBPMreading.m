function [B,T] = SCgetBPMreading(SC,varargin)
% SCgetBPMreading
% ===============
%
% NAME
% ----
% SCgetBPMreading - Calculates BPM readings based on the current injection scheme.
%
% SYNOPSIS
% --------
% `[B, T] = SCgetBPMreading(SC [, options])`
%
%
% DESCRIPTION
% -----------
% This function calculates the particel trajectories based on the current
% injection setup as defined in `SC.INJ` and calculates the corresponding BPM
% readings. The injection setup is a structure with the fields:
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
%	String defining the tracking mode. If set to orbit mode ('ORB'), the AT
%	function `findorbit6` is used to calculate the trajectories. Otherwise,
%	bunches are generated and tracking is performed.  In both cases the
%	corresponding BPM readings are calculated. If the tracking mode is
%	'pORB', the pseudo-orbit is calculated by averaging the turn-by-turn
%	BPM readings.
%
% If the global variable `plotFunctionFlag` is 1, the tracking results are plotted.
%
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
% `'BPMords'` (`SC.ORD.BPM`):: List of BPM ordinates at which the reading should be returned
%
% RETURN VALUES
% -------------
% `B`::
%	BPM readings
% `T`::
%	Particle trajectories
%
% GLOBALS
% -------
% `plotFunctionFlag`:: If true, each BPM reading is plotted
%
% EXAMPLES
% --------
% Switch on plotting, track 100 particles for 2 turns and store the BPM readings in `B`.
% ------------------------------------------------------------------
% global plotFunctionFlag
% plotFunctionFlag  = 1;
% SC.INJ.trackMode  = 'TBT';
% SC.INJ.nParticles = 100;
% SC.INJ.nTurns     = 2;
% B = SCgetBPMreading(SC);
% ------------------------------------------------------------------
%
% SEE ALSO
% --------
% *SCgenBunches*, *SCregisterBPMs*, *SCplotBPMreading*


	global plotFunctionFlag
	
	% Parse optional arguments
	p = inputParser;
	addOptional(p,'BPMords',[]);
	parse(p,varargin{:});

	
	% Switch orbit or tracking mode
	if strcmp(SC.INJ.trackMode,'ORB')
		nTurns     = 1;
		nParticles = 1;
	else
		nTurns     = SC.INJ.nTurns;
		nParticles = SC.INJ.nParticles;
	end

	% Prelocate BPM readings
	B1 = nan(2,nTurns*length(SC.ORD.BPM),SC.INJ.nShots);

	% Check if plotting is desired
	if plotFunctionFlag
		T1          = nan(6,nTurns*nParticles*length(SC.RING),SC.INJ.nShots);
		refOrds     = 1:length(SC.RING);
	else
		refOrds     = SC.ORD.BPM;
	end

	% Loop over injections
	for nShot=1:SC.INJ.nShots

		% Switch orbit or tracking mode
		if strcmp(SC.INJ.trackMode,'ORB')
			% Calculate closed orbit
			T = findorbit6(SC.RING, refOrds);%,SC.INJ.Z0);
		else
			% Generate initial particle distribution
			Zin = SCgenBunches(SC);

			% Calculate trajectories
			T = atpass(SC.RING, Zin, 1, nTurns, refOrds);
		end

		% Account for particle lost (AT sets only x to nan)
		T(:,isnan(T(1,:))) = nan;

		if plotFunctionFlag
			% Store for plotting
			T1(:,:,nShot) = T;
		end

		% Get BPM reading
		B1(:,:,nShot) = calcBPMreading(SC,T,'atAllElements',plotFunctionFlag);
	end

	% Calculate average BPM reading
	B = mean(B1,3,'omitnan');

	% Plot trajectories
	if plotFunctionFlag
		SCplotBPMreading(SC,B,T1);
	end

	% Get pseudo-orbit BPM reading
	if strcmp(SC.INJ.trackMode,'pORB')
		Bpseudo = nan(2,length(SC.ORD.BPM));
		for nBPM=1:length(SC.ORD.BPM)
			Bpseudo(:,nBPM) = mean(B(:,nBPM:length(SC.ORD.BPM):end),2,'omitnan');
		end
		B = Bpseudo;
	end

	% Check if user specified list of BPMs has been defined
	if ~isempty(p.Results.BPMords)
		% Select BPMs for return array
		ind = find(ismember(SC.ORD.BPM,p.Results.BPMords));
		% Check if all ordinates have been recognized
		if length(ind)~=length(p.Results.BPMords)
			warning('Not all specified ordinates are registered BPMs.')
		end
		% Get multiturn indices
		if strcmp(SC.INJ.trackMode,'TBT')
			ind = [0:(nTurns-1)]*length(SC.ORD.BPM)+ind';
		end
		% Define final array of BPM readings
		B = B(:,ind(:));
	end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary function

function B = calcBPMreading(SC,T,varargin)
% Calculates the BPM readings based on a multi-{particle,turn,shot} trajectories. The beam is
% considered lost if the relative amount of lost particles exceeds the value specified in
% `SC.INJ.beamLostAt` and the corresponding BPM reading is `NaN`.
%
% INPUTS
% ------
% `SC`::   SC base structure.
% `T`::    Particle trajectories, e.g. as given by `atpass`.
%
% OPTIONS
% -------
% 'atAllElements' (0)::
%   Flags if input trajectories are given only at BPM ordinates or at all lattice elements.
%
% RETURN VALUES
% -------------
% `B`::
%	[2 x (length(BPMords)*nTurns)] array of x and y BPM readings [m]


	% Parse optional arguments
	p = inputParser;
	addOptional(p,'atAllElements',0);
	parse(p,varargin{:});

	% Switch orbit or tracking mode
	if strcmp(SC.INJ.trackMode,'ORB')
		nTurns     = 1;
		BPMnoise   = cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'NoiseCO'))';
		nParticles = 1;
	else
		nTurns     = SC.INJ.nTurns;
		BPMnoise   = cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'Noise'))';
		nParticles = SC.INJ.nParticles;
	end

	
	
	% Define BPM errors along trajectory (multiple turns)
	BPMoffset   = repmat(cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'Offset'))' + cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'SupportOffset'))',1,nTurns);
	BPMcalError = repmat(cell2mat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'CalError'))',1,nTurns);
	BPMroll     = repmat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'Roll')',1,nTurns) + repmat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'SupportRoll')',1,nTurns);
	BPMnoise    = repmat(    BPMnoise   ,1,nTurns) .* SCrandnc(2,2,nTurns*length(SC.ORD.BPM));
	BPMsumError = repmat(atgetfieldvalues(SC.RING(SC.ORD.BPM),'SumError')',1,nTurns);
	
	% Check if trajectories were given not only at BPMs
	if p.Results.atAllElements
		% Indices of BPMs in trajectories
		nE = reshape((0:nTurns-1)*length(SC.RING)+SC.ORD.BPM',1,[]);
	else
		nE = 1:(length(SC.ORD.BPM)*nTurns);
	end
	
	% Check for multiparticle tracking
	if nParticles > 1

		% Reshape trajectories in tensor
		M = SCparticlesIn3D(T,nParticles);

		% Read x and y positions
		Tx = squeeze(M(1,nE,:));
		Ty = squeeze(M(3,nE,:));

		% Calculate center of charge
		Bx1 = mean(Tx,2,'omitnan')';
		By1 = mean(Ty,2,'omitnan')';
		
		% Check for detactable BPM signal
		beamLost  = find( sum( isnan( Tx ),2 )' .* (1 + BPMsumError.*SCrandnc(2,size(BPMsumError))) > ( nParticles * SC.INJ.beamLostAt ),1);

		% Reflect effective non detectable signal
		Bx1(beamLost:end) = nan;
		By1(beamLost:end) = nan;

	else
		% Read x and y positions
		Bx1 = T(1,nE);
		By1 = T(3,nE);
	end

	% Add roll error
	Bx = cos(BPMroll) .* Bx1 - sin(BPMroll) .* By1;
	By = sin(BPMroll) .* Bx1 + cos(BPMroll) .* By1;

	% Add calibration and offset
	Bx = (Bx - BPMoffset(1,:)) .* (1+BPMcalError(1,:));
	By = (By - BPMoffset(2,:)) .* (1+BPMcalError(2,:));

	% Add noise
	Bx = Bx + BPMnoise(1,:);
	By = By + BPMnoise(2,:);

	% Store BPM reading in output variable
	B(1,:) = Bx;
	B(2,:) = By;

end


