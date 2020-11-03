function [DA,RMAXs,thetas] = SCdynamicAperture(RING,dE,varargin)
% SCdynamicAperture
% =================
%
% NAME
% ----
% SCdynamicAperture - calculates the dynamic aperture of a ring
%
% SYNOPSIS
% --------
% `[DAs, RMAXs, thetas] = SCdynamicAperture(RING, dE [, options])`
%
%
% DESCRIPTION
% -----------
% Calculates the dynamic aperture (i.e. the area of stable particle motion) of
% `RING` at energy `dE`. The general strategy is to find the
% maximum radii at which the particle motion is stable along a number of
% straight lines, leaving the origin at angles `thetas`. The dynamic aperture
% is then approximated as the area of the resulting polygon.
%
%
% INPUT
% -----
% `RING`::
%	Lattice cell structure.
% `dE`::
%	Momentum deviation.
%
% OPTIONS
% -------
% The following options can be specified as name-value pairs:
%
% `'bounds'` (`[0,1e-3]`)::
%	a 1x2 array containing a best guess for the upper and lower boundaries
%	of the maximum radius. These inital boundaries are automatically
%	refined in this routine, so a rough guess is good enough.
% `'nturns'` (1000)::
%	number of turns used to determine whether a particle is stable.
% `'thetas'` (`linspace(0,2*pi,16)`)::
%	angles at which the maximum radii are evaluated.
% `'accuracy'` (`1e-4`)::
%	is the accuracy to which the dynamic aperture is determined.
% `'launchOnOrbit'` (0)::
%	If true, particles are launched on closed orbit (findorbit4), otherwise on axis
% `'centerOnOrbit'` (1)::
%	If true, the closed orbit (findorbit4) is subtracted from the DA coordinates, which 
%   is advised for a corrected machine with support structure misalignments.
% `'auto'` (0)::
%	if >0, this number of automatically determined sampling points is used,
%	taking into account a presumed near-elliptical shape of the DA. In this
%	case `'thetas'` is ignored.
% `'plot'` (0)::
%	if true, progress is plotted.
% `'verbose'` (0)::
%	if true, debug messages are printed.
%
% GLOBALS
% -------
% `runParallel`:: If true, a parfor loop is executed instead of a regular for loop.
%
% RETURN VALUES
% -------------
% `DA`::
%	Dynamic aperture in m^2.
% `RMAXs`::
%	Maximum radii at the evaluated `thetas`. `length(thetas)` array.
% `thetas`::
%	Angles at which the maximum radii were evaluated.


	% Parse input
	p = inputParser;
	addOptional(p,'bounds',[0,1e-3]);
	addOptional(p,'nturns',1000);
	addOptional(p,'nsteps',0);
	addOptional(p,'thetas',linspace(0,2*pi,16));
	addOptional(p,'accuracy',1e-4);
	addOptional(p,'launchOnOrbit',0);
	addOptional(p,'centerOnOrbit',1);
	addOptional(p,'auto',0);
	addOptional(p,'plot',0);
	addOptional(p,'verbose',0);
	parse(p,varargin{:});
	par = p.Results;

	inibounds = par.bounds;
	nturns = par.nturns;
	thetas = par.thetas;
	
	% Check if parallel computation should be used
	global runParallel
	if runParallel
		parforArg = Inf;
	else
		parforArg = 0;
	end

	% Issue warning if nsteps was defined
	if par.nsteps~=0; warning('nsteps no longer supported; continuing with binary search.'); end;

	% If requested, auto generate thetas
	if par.auto>0; [~,thetas] = autothetas(RING,dE,par.auto,varargin{:}); end;


	[~,sidx] = sort(abs(inibounds)); % Sort bounds w.r.t absolute value
	inibounds = inibounds(sidx);

	ZCO = zeros(6,1);
	if par.launchOnOrbit
		tmp = findorbit4(RING,0,[1]); % Closed orbit at reference points
		if ~isnan(tmp(1))
			ZCO(1:4) = tmp;
		end
	end
	ZCO(5) = dE;
	
	RMAXs = nan(length(thetas),1); % Initialize output array
	DA = nan;

	parfor (cntt = 1:length(thetas),parforArg) % Loop over angles
		theta=thetas(cntt);

		bounds = inibounds;

		fatpass(RING,nan(6,1),1,1,[1]); % Fake Track to initialize lattice

		% Scale boundaries up until maxr is included
		scales=0;
		while scales<16
			if check_bounds(RING,ZCO,nturns,theta,bounds); break; end;
			bounds = scale_bounds(bounds,10);
			scales = scales + 1;
			if par.verbose; fprintf('Scaled: %e %e\n',bounds(1),bounds(2)); end;
		end

		% Refine boundaries until requested accuracy is reached
		while abs((bounds(2)-bounds(1))/max(bounds)) > par.accuracy
			bounds = refine_bounds(RING,ZCO,nturns,theta,bounds);
			if par.verbose; fprintf('Refined: %e %e\n',bounds(1),bounds(2)); end;
		end

		RMAXs(cntt)=mean(bounds); % Store mean of final boundaries

	end
	if par.plot
		figure(6232);
		scatter(cos(thetas)'.*RMAXs,sin(thetas)'.*RMAXs);
		set(gca,'xlim',18E-3*[-1 1],'ylim',18E-3*[-1 1])
		drawnow;
	end

	% Calculate DA area
	dthetas = diff(thetas);
	r0 = RMAXs(1:(end-1));
	r1 = RMAXs(2:end);
	DA = sum(sin(dthetas) .* r0' .* r1' / 2.);

	% Center DA around closed orbit
	if par.centerOnOrbit
		tmp = findorbit4(RING,0,1);
		if ~isnan(tmp(1))
			[x,y] = pol2cart(thetas,RMAXs');
			x = x - tmp(1);
			y = y - tmp(3);
			[thetas,RMAXs] = cart2pol(x,y);
			RMAXs = RMAXs';
		end
	end
end


function res = check_bounds(RING,ZCO,nturns,theta,boundsIn)
	rmin = boundsIn(1);
	rmax = boundsIn(2);

	Zmin = ZCO;
	Zmin(1) = rmin * cos(theta);
	Zmin(3) = rmin * sin(theta);
	Zmax = ZCO;
	Zmax(1) = rmax * cos(theta);
	Zmax(3) = rmax * sin(theta);

	ROUT = fatpass(RING,[Zmin,Zmax],0,nturns,[1]); % Track
	RLAST = ROUT(:,end-1:end); % Get positions after last turn
	if (~isnan(RLAST(1,1)) && isnan(RLAST(1,2)))
		res = true;
	else
		res = false;
	end
end


function bounds = refine_bounds(RING,ZCO,nturns,theta,boundsIn)
	rmin = boundsIn(1);
	rmax = boundsIn(2);

	rmid = mean(boundsIn);
	Z = ZCO;
	Z(1) = rmid * cos(theta);
	Z(3) = rmid * sin(theta);

	ROUT = fatpass(RING,Z,0,nturns,[1]); % Track

	RLAST = ROUT(:,end); % Get positions after last turn

	if isnan(RLAST(1)) % Midpoint is outside DA
		bounds = [rmin,rmid];
	else % Midpoint is within DA
		bounds = [rmid,rmax];
	end
end

function out = scale_bounds(bounds,alpha)
	lower = mean(bounds)-(mean(bounds)-bounds(1)) * alpha;
	upper = mean(bounds)-(mean(bounds)-bounds(2)) * alpha;
	% Prohibit zero-crossing during scaling
	if sign(lower) ~= sign(bounds(1))
		lower = 0.0;
	end
	if sign(upper) ~= sign(bounds(2))
		upper = 0.0;
	end
	out = [lower,upper];
end

function r = fatpass(varargin)
	% fatpass - fake atpass to circumvent some madlab parfor shizzle. I don't care anymore.
	r = atpass(varargin{:});
end

function [rout,tout] = autothetas(RING,dE,nt,varargin)
	tin = linspace(0,2*pi*3/4,4);
	[DA,rs,ts] = SCdynamicAperture(RING,dE,varargin{:},'thetas',tin,'auto',0);
	a=(rs(1)+rs(3))/2;
	b=(rs(2)+rs(4))/2;
	e=sqrt(1-b^2/a^2);
	mu=acosh(1/e);
	nu=linspace(0,2*pi,nt);
	rout = abs(cosh(mu+1i*nu));
	tout = angle(cosh(mu+1i*nu));
end
