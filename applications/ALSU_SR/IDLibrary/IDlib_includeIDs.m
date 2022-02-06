function RING = IDlib_includeIDs(RING,IDs,varargin)
% IDlib_includeIDs
% ================
%
% NAME
% ----
% IDlib_includeIDs - Includes IDs in ALSU-SR plain lattice
%
% SYNOPSIS
% --------
% `RING = IDlib_includeIDs(RING,IDs [,options])`
%
% DESCRIPTION
% -----------
% This function includes a list of IDs at various center of straights into the provided lattice.
% Possible IDs are 'LEDA','COSMIC','TENDER','EPU36','XType','EPU50','EPU35','EPU90' (all as kick
% maps) and 'U114' (series of SBENDs).
%
% INPUT
% -----
% `RING`::
% 	AT lattice cell structure structure.
% `IDs`::
% 	Cell string with ID names.
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'Section'` ([])::
%	Section of ID.
% `factor` (`ones(size(IDs))`)::
%	Scaling factor if ID (only for series of SBENDs)
%
% RETURN VALUE
% ------------
% `RING`::
% 	AT lattice cell structure structure with IDs.
%
% SEE ALSO
% --------
% *IDLib_closeID*

	p = inputParser;
	addOptional(p,'Section',[]);
	addOptional(p,'factor',ones(size(IDs)));
	parse(p,varargin{:});
	par = p.Results;	
	
	for nID=1:length(IDs)
		% Load ID elements
		IDelements = getIDElements(RING,IDs{nID},par.factor(nID));
		
		% Include ID elements
		RING = includeIDElements(RING,IDelements,par.Section(nID));
	end
	
end

function IDelements = getIDElements(RING,name,factor)
	energy = 2E9;
	
	switch name
		case {'LEDA','COSMIC','TENDER','EPU36','XType','EPU38','EPU50','EPU35','EPU90'}
			% Create placeholder for kick map IDs 
			IDelements.FamName    = name;
			IDelements.Length     = 0;
			IDelements.PassMethod = 'IdentityPass';
			IDelements.Energy     = energy;
			
		case {'U114'}
			% Create series of SBENDs to model wiggler
			global FAMLIST 
			FAMLIST = cell(0);
			
			tmp = SCgetOrds(RING,'D11$');
			dw11 = RING{tmp(1)}.Length;
			
			% Define start and ned drift and shikane elements
			DW11      = drift('DWr11',   (8*dw11-29*0.114)/2,  'DriftPass');
			W11N1     = @(name,th) sbend(name, 0.114/4, th, 0, th, 0,  'DriftPass');
			W11N2     = @(name,th) sbend(name, 0.114/4, th, th, 0, 0,  'DriftPass');
			W11S1     = @(name,th) sbend(name, 0.114/4, -th, -th, 0, 0,  'DriftPass');
			W11S2     = @(name,th) sbend(name, 0.114/4, -th, 0, -th, 0,  'DriftPass');
			W11period = @(name,th) [W11N1(name,th), W11S1(name,th), W11S2(name,th), W11N2(name,th)];
			
			% Wiggler angle
			WigThetaU114 = sqrt(sum((1.9.*sin(0:0.001:2*pi)).^2)./length(0:0.001:2*pi))*(0.114./4)/(2.0/0.2998);
	
			% Define total wiggler
			wiggler = [DW11 repmat(W11period(name,WigThetaU114/factor),1,29) DW11];
			
			% Generate lattice elements
			IDelements=cell(size(wiggler));
			for i=1:length(wiggler)
				IDelements{i} = FAMLIST{wiggler(i)}.ElemData;
				IDelements{i}.Class  = 'Drift';
				IDelements{i}.Energy = energy;
			end

		otherwise
			error(fprintf('Unsupported ID name: ''%s''.',name))
	end
end



function RING = includeIDElements(RING,IDelements,section)
	
	% Split RING to include IDs
	CoSord = SCgetOrds(RING,sprintf('CenOfStr%d$',section));
	
	if section==12
		if length(IDelements)==1
			warning('Not implemented yet.')
		else
			warning('Works only if there are 8 D11 elements per straight!')

			RING = [IDelements(length(IDelements)/2+1:end) RING(6:end-4) IDelements(1:length(IDelements)/2) ];
		end
	else
		if length(IDelements)==1
			RINGup   = RING(1:CoSord);
			RINGdown = RING((CoSord+1):end);
		else
			driftOrds = SCgetOrds(RING,'D11$');
			
			cutOut = [driftOrds(find(driftOrds<CoSord,4,'last')) driftOrds(find(driftOrds>CoSord,4))];
			
			RINGup   = RING(1:cutOut(1)-1);
			RINGdown = RING((cutOut(end)+1):end);
		end
		
		RING = [RINGup IDelements RINGdown];
	end
end
