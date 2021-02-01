function res = SCgetSupportOffset(SC,s)
% SCgetSupportOffset
% ==================
%
% NAME
% ----
% SCgetSupportOffset - Calculates the combined support structure offset at a certain location
%
% SYNOPSIS
% --------
% `res = SCgetSupportOffset(SC, s)`
%
%
% DESCRIPTION
% -----------
% This function evaluates the total offset of the support structures that have
% been defined via *SCregisterSupport* at the longitudinal positions `s` by linearly interpolating 
% between girder start- and endpoints (which may be mounted on sections and plinths).
% Note that hereby longitudinal offsets of support structures are only correctly calculated in the 
% case of straight girders (no magnet with bending angle on it) and with equal longitudinal offsets
% at their start- and endpoints.
%
% INPUTS
% ------
% `SC`::
%	The `SC` core structure.
% `s`::
%	Array of s-positions at which the offset is evaluated.
%
% RETURN VALUE
% ------------
% `res`::
%	`[3,length(s)]`-array containing the x/y/z offsets at `s`.
%
% SEE ALSO
% --------
% *SCregisterSupport*, *SCplotSupport*

	% TODO: inlcude longitudinal support structure offsets

	% No girders -> no cookies
	if ~isfield(SC.ORD,'Girder')
		res=zeros(2,length(s));
		return;
	end

	% Read elemet lengths from RING
	lengths = nan(1,length(SC.RING));
	for ord=1:length(SC.RING)
		if isfield(SC.RING{ord},'Length')
			lengths(ord)=SC.RING{ord}.Length;
		else
			lengths(ord)=0.0;
		end
	end


	C=sum(lengths); % Circumference
	N=length(SC.RING); % Number of elements
	sposEND=cumsum(lengths); % s-position of element end points
	sposMID=sposEND-lengths./2; % s-sposition of element middles


	% Generate arrays containing the x/y/z-offsets of
	% the beginning/end of the girders
	aG=SC.ORD.Girder(1,:); % Beginning ordinates
	bG=SC.ORD.Girder(2,:); % End ordinates
	sGa=sposMID(aG);
	sGb=sposMID(bG);
	for i=1:length(aG) % Read girder offsets from elements
		xGa(i)=SC.RING{aG(i)}.GirderOffset(1);
		yGa(i)=SC.RING{aG(i)}.GirderOffset(2);
		xGb(i)=SC.RING{bG(i)}.GirderOffset(1);
		yGb(i)=SC.RING{bG(i)}.GirderOffset(2);
		zGa(i)=SC.RING{aG(i)}.GirderOffset(3);
		zGb(i)=SC.RING{bG(i)}.GirderOffset(3);
	end


	% Add the Plinth-offsets to the girder offsets
	if isfield(SC.ORD,'Plinth')
		aP=SC.ORD.Plinth(1,:); % Beginning ordinates
		bP=SC.ORD.Plinth(2,:); % End ordinates

	
		% Plinth s-positions
		sPa=sposMID(aP);
		sPb=sposMID(bP);

		% A/B are matrices:
		% 1: if mod(aG_i,N) \in [aP_j,bP_j]
		% 0: else
		A=mod_int_cont(aG,N,aP,bP);
		B=mod_int_cont(bG,N,aP,bP);
		
		% Read section offsets from elements
		for i=1:length(aP) 
			xPa(i) = SC.RING{aP(i)}.PlinthOffset(1);
			yPa(i) = SC.RING{aP(i)}.PlinthOffset(2);
			zPa(i) = SC.RING{aP(i)}.PlinthOffset(3);
			xPb(i) = SC.RING{bP(i)}.PlinthOffset(1);
			yPb(i) = SC.RING{bP(i)}.PlinthOffset(2);
			zPb(i) = SC.RING{bP(i)}.PlinthOffset(3);
		end
			
		% Interpolate between plinth start- and endpoints and add offset to girders which are on plinths
		xGa = xGa + limp(sGa,C,sPa,xPa,sPb,xPb)';
		xGb = xGb + limp(sGb,C,sPa,xPa,sPb,xPb)';
		yGa = yGa + limp(sGa,C,sPa,yPa,sPb,yPb)';
		yGb = yGb + limp(sGb,C,sPa,yPa,sPb,yPb)';
		zGa = zGa + limp(sGa,C,sPa,zPa,sPb,zPb)';
		zGb = zGb + limp(sGb,C,sPa,zPa,sPb,zPb)';
	end
	


	% Add the section-offsets to the girder(+plinth) offsets
	if isfield(SC.ORD,'Section')
		aS=SC.ORD.Section(1,:); % Beginning ordinates
		bS=SC.ORD.Section(2,:); % End ordinates
		% Section s-positions
		sSa=sposMID(aS);
		sSb=sposMID(bS);

		% Read section offsets from elements
		for i=1:length(aS) 
			xSa(i) = SC.RING{aS(i)}.SectionOffset(1);
			ySa(i) = SC.RING{aS(i)}.SectionOffset(2);
			zSa(i) = SC.RING{aS(i)}.SectionOffset(3);
			xSb(i) = SC.RING{bS(i)}.SectionOffset(1);
			ySb(i) = SC.RING{bS(i)}.SectionOffset(2);
			zSb(i) = SC.RING{bS(i)}.SectionOffset(3);
		end
			
		% Interpolate between section start- and endpoints and add offset of sections to girders (+plinths)
		xGa = xGa + limp(sGa,C,sSa,xSa,sSb,xSb)';
		xGb = xGb + limp(sGb,C,sSa,xSa,sSb,xSb)';
		yGa = yGa + limp(sGa,C,sSa,ySa,sSb,ySb)';
		yGb = yGb + limp(sGb,C,sSa,ySa,sSb,ySb)';
		zGa = zGa + limp(sGa,C,sSa,zSa,sSb,zSb)';
		zGb = zGb + limp(sGb,C,sSa,zSa,sSb,zSb)';
	end


	% Calculate final offset at s, by interpolating
	% over the total girder offsets
	dx = limp(s,C,sGa,xGa,sGb,xGb)';
	dy = limp(s,C,sGa,yGa,sGb,yGb)';
	dz = limp(s,C,sGa,zGa,sGb,zGb)';
	res = [dx;dy;dz];
	
end


function res = mod_int_cont(x,C,a,b)
	x=x(:);
	a=a(:)';
	b=b(:)';
	x=mod(x-a,C); % 1-based indexing
	b=mod(b-a,C);
	res = x<b;
end

% function res = limp(x,C,a1,f1,a2,f2)
	% _L_inear _I_nterpolation in _M_odulo space / _P_iecewise
% 	x =x(:);
% 	a1=a1(:)';
% 	f1=f1(:)';
% 	a2=a2(:)';
% 	f2=f2(:)';
% 	x =mod( x-a1,C);
% 	a2=mod(a2-a1,C);
% 
% 	res = f1+(f2-f1)./a2.*x;
% 	res = res .* (x<a2);
% end

function res = limp(x,C,a1,f1,a2,f2)
	% _L_inear _I_nterpolation in _M_odulo space / _P_iecewise
	x =x(:)';
	a1=a1(:)';
	f1=f1(:)';
	a2=a2(:)';
	f2=f2(:)';
	res = zeros(length(x),1);
	
	for n=1:length(a1)
		if a1(n)<a2(n)
			ind = intersect(find(x>=a1(n)),find(x<=a2(n)));
			res(ind) = interp1([a1(n) a2(n)],[f1(n) f2(n)],x(ind));
		else
			ind1 = find(x<=a2(n));
			ind2 = find(x>=a1(n));
			res(ind1) = interp1([a1(n) a2(n)+C],[f1(n) f2(n)],C+x(ind1));
			res(ind2) = interp1([a1(n) a2(n)+C],[f1(n) f2(n)],x(ind2));
		end
	end

end
