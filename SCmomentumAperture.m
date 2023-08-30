function dbounds = SCmomentumAperture(RING,REFPTS,inibounds,varargin)
% SCmomentumAperture
% ==================
%
% NAME
% ----
% SCmomentumAperture - Calculates the momentum aperture
%
% SYNOPSIS
% --------
% `dbounds = SCmomentumAperture(RING, REFPTS, inibound [, options])`
%
%
% DESCRIPTION
% -----------
% Calculates the momentum aperture of `RING` at the given `REFPTS`.  `dbounds`
% is a `[2,length(REFPTS)]`-array containing the upper and lower bounds of the
% local momentum aperture either the positive (default) or negative direction,
% depending on the sign of `'inibounds'`, which is an initial guess for the
% bounds of the momentum apeture. This is automatically refined in this
% routine, so a rough guess is good enough.
%
%
% INPUT
% -----
% `RING`::
%	Lattice cell structure.
% `REFPTS`::
%	`[1xN]` array of lattice ordinates at which the MA should be calculated.
% `inibounds`::
%	`[inner,outer]` array of inner and outer bounds intial guess for the MA.
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'nturns'` (1000)::
%	is the number of turns used to determine whether a particle is stable.
% `'accuracy'`(1e-4)::
%	is the accuracy to which the momentum aperture is determined.
% `'stepsize'`(1e-3)::
%	stepsize, with which the boundaries are expanded, when the critial
%	point is not within them.
% `'debug'` (0)::
%	if true, debug information is printed.
% `'plot'` (0)::
%	if true, results are plotted.
%
% GLOBALS
% -------
% `runParallel`:: If true, a parfor loop is executed instead of a regular for loop.
%
% RETURN VALUE
% ------------
% `dbounds`:: `2xlength(REFPTS)` array containing the local momentum apertures
%	 bounds.
%
%
% EXAMPLES
% --------
% Calculates the positive MA bounds for the first 10 lattice elements of `SC.RING` with 
% intial guess of [0,4E-2] over 1000 turns. Debug information is printed and the results
% are plotted.
% ------------------------------------------------------------------
% dbounds = SCmomentumAperture(SC.RING,1:10,[0,4E-2],'debug',1,'plot',1);
% ------------------------------------------------------------------
%
% Calculates the negative MA bounds at all `SD` magnets of `SC.RING` with the initial
% guess of [-5E-3,-2E-2] over 10000 turns in parallel mode. 
% ------------------------------------------------------------------
% ords = SCgetOrds(SC.RING,'SD');
% global runParallel
% runParallel = 1;
% dbounds = SCmomentumAperture(SC.RING,ords,[-5E-3,-2E-2],'nturns',10000);
% ------------------------------------------------------------------



	% Parse input
	p = inputParser;
	addOptional(p,'nturns',1000);
	addOptional(p,'accuracy',1e-4);
	addOptional(p,'stepsize',1e-3);
	addOptional(p,'plot',0);
	addOptional(p,'debug',0);
	parse(p,varargin{:});
	par = p.Results;

	
	% Check if parallel computation should be used
	global runParallel
	if runParallel
		parforArg = Inf;
	else
		parforArg = 0;
	end

	
	% Initialize output arrays
	dboundHI = zeros(length(REFPTS),1); 
	dboundLO = zeros(length(REFPTS),1); 
	
	% Closed orbits at reference points
	ZCOs = findorbit6(RING,REFPTS); 
	if any(~isfinite(ZCOs(:)))
		dbounds = [dboundHI,dboundLO];
		warning('Closed Orbit could not be determined during MA evaluation. MA set to zero.');
		return;
	end

	% Loop over reference points
	parfor (i=1:length(REFPTS),parforArg) 
		ord = REFPTS(i);

		% Reset bounds to initial values
		local_bounds=inibounds;

		% Shift so that ord-th element is the first
		SHIFTRING = circshift(RING,-ord+1);

		% Closed orbit at current element
		ZCO = ZCOs(:,i); 

		% Scale boundaries up until dp is included
		while ~check_bounds(local_bounds,SHIFTRING,ZCO,par.nturns)
			local_bounds = increment_bounds(local_bounds,par.stepsize);
			if par.debug; fprintf('ord: %d; Incremented: %+0.5e %+0.5e\n',ord,local_bounds(1),local_bounds(2)); end
		end

		% Refine boundaries until requested accuracy is reached
		while abs((local_bounds(2)-local_bounds(1))/max(local_bounds)) > par.accuracy
			local_bounds = refine_bounds(local_bounds,SHIFTRING,ZCO,par.nturns);
			if par.debug; fprintf('ord: %d; Refined: %e %e\n',ord,local_bounds(1),local_bounds(2)); end
		end

		% Store final boundaries
		dboundHI(i) = local_bounds(1);
		dboundLO(i) = local_bounds(2);

		if par.debug; fprintf('ord: %d; Found: %+0.5e %+0.5e\n',ord,local_bounds(1),local_bounds(2)); end

	end
	
	% Store function output
	dbounds = [dboundHI,dboundLO];
		
	% Plot results
	if par.plot
		spos=findspos(RING,REFPTS);
		figure(81222);clf;hold on;
		plot(spos,dboundHI,'kx-');
		plot(spos,dboundLO,'rx-');
		xlabel('s [m]');ylabel('MA');drawnow;
	end

end

function local_bounds = refine_bounds(local_bounds,RING,ZCO,nturns)
	% bounds shall be absolute-ordered i.e.
	% [0.1,0.2] or [-0.1,-0.3]
	dmean = mean(local_bounds);

	% Add momentum deviations to closed orbit
	Z0 = ZCO;
	Z0(5) = Z0(5) + dmean;

	% Track
	ROUT = atpass(RING,Z0,1,nturns,[]);

	if isnan(ROUT(1)) % Particle dead :(
		local_bounds(2) = dmean; % Set abs-upper bound to midpoint
	else % Particle alive :)
		local_bounds(1) = dmean; % Set abs-lower bound to midpoint
	end

end

function check = check_bounds(local_bounds,RING,ZCO,nturns)
	% Returns true when the lower (upper) bound is (un)stable

	% Add momentum deviations to closed orbit
	Z=[ZCO,ZCO];
	Z(5,:)=Z(5,:)+local_bounds(:)';

	% Track
	ROUT = atpass(RING,Z,1,nturns,[]);

	% Throw error if the less-deviating particle is unstable and the other is not
	if isnan(ROUT(1,1)) && ~isnan(ROUT(1,2))
		warning('Closer-to-momentum particle is unstable. This shouldnt be!');
	end

	% Closer-to-momentum should be stable,
	% farther-from-momentum particle should be unstable
	check = ~isnan(ROUT(1,1)) && isnan(ROUT(1,2));

end

function local_bounds = increment_bounds(local_bounds,stepsize)
	local_bounds = local_bounds + [-1,1] .* mysign(local_bounds)*stepsize;
end

function s = mysign(v)
	s = 1 - 2 * (v<0);
end

function out = scale_bounds(local_bounds,alpha)
	lower = mean(local_bounds)-(mean(local_bounds)-local_bounds(1)) * alpha;
	upper = mean(local_bounds)-(mean(local_bounds)-local_bounds(2)) * alpha;
	% Prohibit zero-crossing during scaling
	if sign(lower) ~= sign(local_bounds(1))
		lower = 0.0;
	end
	if sign(upper) ~= sign(local_bounds(2))
		upper = 0.0;
	end
	out = [lower,upper];
end
