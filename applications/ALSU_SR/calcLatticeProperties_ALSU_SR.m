function OUT = calcLatticeProperties_ALSU_SR(SC,varargin)
% calcLatticeProperties_ALSU_SR
% =============================
%
% NAME
% ----
% calcLatticeProperties_ALSU_SR - Calculates various lattice properties
%
% SYNOPSIS
% --------
% `OUT = calcLatticeProperties_ALSU_SR(SC [, options])`
%
%
% DESCRIPTION
% -----------
% This function compares various properties between the lattices `SC.RING` and
% `SC.IDEALRING`.
%
%
% RETURN VALUE
% ------------
% `OUT` is a structure with the fields
%
% `orbit`::
%	RMS orbit difference in horizontal/vertical direction
% `betaBeat`::
%	RMS beta difference in horizontal/vertical direction
% `dispBeat`::
%	RMS dispersion difference in horizontal/vertical direction
% `tuneShift`::
%	Tune difference in horizontal/vertical direction
% `Chromaticity`::
%	Chomaticity of `SC.RING`
% `Coupling`::
%	Frobenius norms of the off-diagonal one-turn 4x4 transfer matrix of
%	`SC.RING`
% `Emittance`::
%	Emittance of `SC.RING` in horizontal/vertical direction
% `DA`, `RMAXs`, `thetas`::
%	Dynamic aperture `SC.RING` and the angles and maximum
%	radii used for the evalutation


	p = inputParser;
	addOptional(p,'DAsteps',[]);
	addOptional(p,'DAmode','thetas');
	addOptional(p,'sepTunesWithOrds',SCgetOrds(SC.RING,{'^QF1','^QD1'}));
	addOptional(p,'sepTunesDeltaK',[-0.01 0.002]);
	parse(p,varargin{:});
	par = p.Results;

	
	dP = 0;%1e-3;
	
	% Define reference points
	REFPTS = 1:length(SC.RING);
	nREF = length(REFPTS);

	% Calculate RMS orbit
	X0=findorbit6(SC.RING,REFPTS);
	OUT.orbit = sqrt(mean(X0([1 3],:).^2,2));

	% Switch in 4D mode
	RING  = SCcronoff(SC.RING,     'radiationoff','cavityoff'); %
	RING0 = SCcronoff(SC.IDEALRING,'radiationoff','cavityoff'); %

	% Use linopt to get lattice info.
	[ld, nus, xis]  = atlinopt(RING, dP,REFPTS);
	[ld0,nus0, ~ ] = atlinopt(RING0,dP,REFPTS);

	% Generate arrays from arkward linopt output
	beta  = reshape([ld.beta ],[2,nREF]);
	beta0 = reshape([ld0.beta],[2,nREF]);
	dispersion  = reshape([ld.Dispersion ],4,nREF);
	dispersion0 = reshape([ld0.Dispersion],4,nREF);

	dispersion  = dispersion( [1 3],:);
	dispersion0 = dispersion0([1 3],:);

	% Write beta/disp beating (rms of difference)
	OUT.betaBeat = sqrt(mean(((beta-beta0)./beta0).^2,2));
	OUT.dispBeat = sqrt(mean((dispersion-dispersion0).^2,2));

	% Write tune shift (difference)
	OUT.tuneShift = [nus-nus0]';

	% Write chromaticities
	OUT.Chromaticity = xis';


	% Calculate measure of coupling
	% as Frobenius norm of off-diagonal one-turn 4x4 transfer matrix
	% at the injection point.
	M44 = findm44(RING,0,length(RING)+1);
	UPR = M44([1 2],[3 4]);
	LWR = M44([3 4],[1 2]);
	OUT.Coupling(1,1) = sqrt(trace(UPR'*UPR));
	OUT.Coupling(2,1) = sqrt(trace(LWR'*LWR));


	% Dynamic aperture
	if ~isempty(par.DAsteps)
		fprintf('Calculating DA with %d steps\n',par.DAsteps)
		RING  = SCcronoff(SC.RING,'radiationoff','cavityoff');
		[OUT.DA,OUT.RMAXs,OUT.thetas] = SCdynamicAperture(RING,0,'accuracy',1e-6,'plot',0,par.DAmode,linspace(0,2*pi,par.DAsteps));
	end

	
	
	% Emittance % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
	
	% Seperate tunes for fitting
	if ~isempty(par.sepTunesWithOrds) && ~isempty(par.sepTunesDeltaK)
		fprintf('Separating tunes for emittance calculation by changing quadrupoles by %.1e/%.1e\n',par.sepTunesDeltaK)
		for nFam=1:length(par.sepTunesWithOrds)
			SC = SCsetMags2SetPoints(SC,par.sepTunesWithOrds{nFam},2,2,par.sepTunesDeltaK(nFam),'method','add');
		end
	end
	
	
	RING  = SCcronoff(SC.RING,'radiationon','cavityon'); %

	BENDINDEX = findcells(RING,'PassMethod','BndMPoleSymplectic4RadPass');
	QUADSEXTINDEX = findcells(RING,'PassMethod','StrMPoleSymplectic4RadPass');
	RADELEMINDEX = sort([BENDINDEX QUADSEXTINDEX]);

	try
		[ENV, DP, DL] = ohmienvelope(RING(:),RADELEMINDEX, REFPTS);

		sigmas = cat(2,ENV.Sigma);

		RING  = SCcronoff(SC.RING,'radiationoff','cavityoff'); %

		dP = 0.0001;
		eta = 5000*(findorbit4(RING,dP,REFPTS)-findorbit4(RING,-dP,REFPTS));

		ex = ((sigmas(1,:).^2) - (eta(1,:)*DP).^2)./beta(1,:);
		ey = ((sigmas(2,:).^2) - (eta(3,:)*DP).^2)./beta(2,:);

		ey(abs(ey)>10*mean(abs(ey)))=nan;
		ex(abs(ex)>10*mean(abs(ex)))=nan;

		% Write emittance
		OUT.Emittance = [mean(ex,'omitnan');mean(ey,'omitnan')];
	catch

		OUT.Emittance = nan(2,1);
	end
	% E N D Emittance % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % %

	

end






