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
% `SC = SCgetSupportOffset(SC, s)`
%
%
% DESCRIPTION
% -----------
% This function evaluates the total offset of the support structures that have
% been defined via *SCregisterSupport* at the longitudinal positions `s`.
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
%	`[2,length(s)]`-array containing the x/y offsets at `s`.
%
% SEE ALSO
% --------
% *SCregisterSupport*, *SCplotSupport*

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


	% Generate arrays containing the x/y-offsets of
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
	end


	% Add the Plinth-offsets to the girder offsets
	if isfield(SC.ORD,'Plinth')
		aP=SC.ORD.Plinth(1,:); % Beginning ordinates
		bP=SC.ORD.Plinth(2,:); % End ordinates

		numP=length(aP);
		xP=zeros(numP,1);
		yP=zeros(numP,1);
		for i=1:numP % Read plinth offsets from elements
			xP(i) = SC.RING{aP(i)}.PlinthOffset(1);
			yP(i) = SC.RING{aP(i)}.PlinthOffset(2);
		end

		% A/B are matrices:
		% 1: if mod(aG_i,N) \in [aP_j,bP_j]
		% 0: else
		A=mod_int_cont(aG,N,aP,bP);
		B=mod_int_cont(bG,N,aP,bP);

		% Add offset of plinths to girders
		xGa = xGa + (A*xP)';
		xGb = xGb + (B*xP)';
		yGa = yGa + (A*yP)';
		yGb = yGb + (B*yP)';
	end


	% Add the section-offsets to the girder+plinth offsets
	if isfield(SC.ORD,'Section')
		aS=SC.ORD.Section(1,:); % Beginning ordinates
		bS=SC.ORD.Section(2,:); % End ordinates

		numS=length(aS);
		xS=zeros(numS,1);
		yS=zeros(numS,1);
		for i=1:numS % Read section offsets from elements
			xS(i) = SC.RING{aS(i)}.SectionOffset(1);
			yS(i) = SC.RING{aS(i)}.SectionOffset(2);
		end

		% A/B are matrices:
		% 1: if mod(aG_i,N) \in [aS_j,bS_j]
		% 0: else
		A=mod_int_cont(aG,N,aS,bS);
		B=mod_int_cont(bG,N,aS,bS);

		% Add offset of sections to girders+plinths
		xGa = xGa + (A*xS)';
		xGb = xGb + (B*xS)';
		yGa = yGa + (A*yS)';
		yGb = yGb + (B*yS)';
	end


	% Calculate final offset at s, by interpolating
	% over the total girder offsets
	dx = sum(limp(s,C,sGa,xGa,sGb,xGb),2)';
	dy = sum(limp(s,C,sGa,yGa,sGb,yGb),2)';
	res = [dx;dy];

end


function res = mod_int_cont(x,C,a,b)
	x=x(:);
	a=a(:)';
	b=b(:)';
	x=mod(x-a,C); % Sucky 1-based indexing
	b=mod(b-a,C);
	res = x<b;
end

function res = limp(x,C,a1,f1,a2,f2)
	% _L_inear _I_nterpolation in _M_odulo space / _P_iecewise
	x =x(:);
	a1=a1(:)';
	f1=f1(:)';
	a2=a2(:)';
	f2=f2(:)';
	x =mod( x-a1,C);
	a2=mod(a2-a1,C);

	res = f1+(f2-f1)./a2.*x;
	res = res .* (x<a2);
end
