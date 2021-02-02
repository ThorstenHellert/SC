function off = SCgetSupportOffset(SC,s)
% SCgetSupportOffset
% ==================
%
% NAME
% ----
% SCgetSupportOffset - Calculates the combined support structure offset at a certain location
%
% SYNOPSIS
% --------
% `off = SCgetSupportOffset(SC, s)`
%
%
% DESCRIPTION
% -----------
% This function evaluates the total offset of the support structures that have
% been defined via *SCregisterSupport* at the longitudinal positions `s` by linearly interpolating 
% between support structure start- and endpoints (girder + sections + plinths, if registered).
% Note that this calculation may not provide the proper values if magnets with non-zero bending
% angle are within the support structure because it does not account for the rotation of the
% local coordinate system along the beam trajectory.
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
% `off`::
%	`[3,length(s)]`-array containing the [dx/dy/dz] total support structure offsets at `s`.
%
% SEE ALSO
% --------
% *SCregisterSupport*, *SCupdateSupport*, *SCgetSupportRoll*, *SCplotSupport*

	% Initialize output
	off = zeros(3,length(s));
	
	% Read elemet lengths from RING
	lengths = getcellstruct(SC.RING,'Length',1:length(SC.RING));

	% Circumference
	C = sum(lengths);    
	% s-sposition of element center (needed if support structure is defined at element with non-zero length)
	sposMID = cumsum(lengths)-lengths./2; 

	% Loop over support structure types
	for type = {'Section','Plinth','Girder'}
		% Check if support structure is registered
		if isfield(SC.ORD,type{1})
			% Ordinates
			ord1=SC.ORD.(type{1})(1,:); % Beginning ordinates
			ord2=SC.ORD.(type{1})(2,:); % End ordinates
			
			% s-positions
			s1=sposMID(ord1);
			s2=sposMID(ord2);
			
			% Read support structure offsets from elements
			tmpoff1=zeros(3,length(ord1));tmpoff2=zeros(3,length(ord2));
			for i=1:length(ord1)
				tmpoff1(:,i) = SC.RING{ord1(i)}.([type{1} 'Offset']);
				tmpoff2(:,i) = SC.RING{ord2(i)}.([type{1} 'Offset']);
			end
			
			% Interpolate between support structure start- and endpoints and add offset to elements which are on structure
			off(1,:) = off(1,:) + limp(s,C,s1,tmpoff1(1,:),s2,tmpoff2(1,:))';
			off(2,:) = off(2,:) + limp(s,C,s1,tmpoff1(2,:),s2,tmpoff2(2,:))';
			off(3,:) = off(3,:) + limp(s,C,s1,tmpoff1(3,:),s2,tmpoff2(3,:))';
		end
	end
end


function out = limp(x,C,a1,f1,a2,f2)
	% _L_inear _I_nterpolation in _M_odulo space / _P_iecewise

	% Uniform input
	x  = x(:)'; % Evaluation points
	a1 = a1(:)'; % 1st sampling point
	f1 = f1(:)'; % 1st function value
	a2 = a2(:)'; % 2nd sampling point
	f2 = f2(:)'; % 2nd function value
	% Initialize output
	out = zeros(length(x),1);
	
	% Loop over evaluation points
	for n=1:length(a1)
		
		% Check if injection is within sampling points
		if a1(n)<a2(n)
			ind = intersect(find(x>=a1(n)),find(x<=a2(n)));
			out(ind) = interp1([a1(n) a2(n)],[f1(n) f2(n)],x(ind));
		else
			ind1 = find(x<=a2(n));
			ind2 = find(x>=a1(n));
			out(ind1) = interp1([a1(n) a2(n)+C],[f1(n) f2(n)],C+x(ind1));
			out(ind2) = interp1([a1(n) a2(n)+C],[f1(n) f2(n)],x(ind2));
		end
	end
end
