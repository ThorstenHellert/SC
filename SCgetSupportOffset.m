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
	

	% s-positions along the ring
	s0 = findspos(SC.RING,1:length(SC.RING));

	% Initialize offset along the ring
	off0 = zeros(3,length(s0));
	
	% Read elemet lengths from RING
	for n=1:length(SC.RING)
		lengths(n) = SC.RING{n}.Length;
	end
	
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
				tmpoff1(:,i) = off0(:,ord1(i))' + SC.RING{ord1(i)}.([type{1} 'Offset']);
				tmpoff2(:,i) = off0(:,ord2(i))' + SC.RING{ord2(i)}.([type{1} 'Offset']);
			end
			
			% Interpolate between support structure start- and endpoints and add offset to elements which are on structure
			off0(1,:) =  limp(off0(1,:),s0,C,s1,ord1,tmpoff1(1,:),s2,ord2,tmpoff2(1,:))';
			off0(2,:) =  limp(off0(2,:),s0,C,s1,ord1,tmpoff1(2,:),s2,ord2,tmpoff2(2,:))';
			off0(3,:) =  limp(off0(3,:),s0,C,s1,ord1,tmpoff1(3,:),s2,ord2,tmpoff2(3,:))';

		end
	end

	% Workaround for bug in plotting which sometimes occurs when many adjacent markers in lattice file with variing offset
	if ~isequal(s,s0)
		% Get offsets at requested s-positions
		[~,b] = unique(s0);
		
		% Interpolate offset
		off(1,:) = interp1(s0(b),off0(1,b),s,'linear','extrap');
		off(2,:) = interp1(s0(b),off0(2,b),s,'linear','extrap');
		off(3,:) = interp1(s0(b),off0(3,b),s,'linear','extrap');
	else
		% Take offset as piecewise interpolated before
		off = off0;
	end
	
end




function off = limp(off,s,C,s1,ord1,f1,s2,ord2,f2)
	% _L_inear _I_nterpolation in _M_odulo space / _P_iecewise
	
	% off:  current offsets along the ring
	% s:    s-positions along the ring
	% C:    ring circumference
	% s1:   s-position of structure start point
	% ord1: ordinate of structure start point
	% f1:   offset of structure start point
	% s2:   s-position of structure end point
	% ord2: ordinate of structure end point
	% f2:   offset of structure end point
	
	
	% Loop over evaluation points
	for n=1:length(s1)
		if s1(n)==s2(n) % Sampling points have same s-position
			if f1(n)~=f2(n)
				error('Something went wrong.')
			end			
			ind = ord1(n):ord2(n);
			off(ind) = f1(n);			
		elseif s1(n)<s2(n) % Standard interpolation
			ind = ord1(n):ord2(n);
			off(ind) = interp1([s1(n) s2(n)],[f1(n) f2(n)],s(ind),'linear','extrap');
		else % Injection is within sampling points
			ind1 = 1:ord2(n);
			ind2 = ord1(n):length(off);
			
			off(ind1) = interp1([s1(n) s2(n)+C],[f1(n) f2(n)],C+s(ind1),'linear','extrap');
			off(ind2) = interp1([s1(n) s2(n)+C],[f1(n) f2(n)],s(ind2),'linear','extrap');
		end
	end
end



% function out = limp(x,C,a1,f1,a2,f2)
% 	% _L_inear _I_nterpolation in _M_odulo space / _P_iecewise
% 
% 	% Uniform input
% 	x  = x(:)'; % Evaluation points
% 	a1 = a1(:)'; % 1st sampling point
% 	f1 = f1(:)'; % 1st function value
% 	a2 = a2(:)'; % 2nd sampling point
% 	f2 = f2(:)'; % 2nd function value
% 	% Initialize output
% 	out = zeros(length(x),1);
% 	
% 	% Loop over evaluation points
% 	for n=1:length(a1)
% 		
% 		
% 		if a1(n)==a2(n) % Sampling points have same s-position
% 			if f1~=f2
% 				error('Something went wrong.')
% 			end
% 			ind = find(x==a1(n));
% 			out(ind) = f1(n);
% 		elseif a1(n)<a2(n) % Standard interpolation
% 			ind = intersect(find(x>=a1(n)),find(x<=a2(n)));
% 			out(ind) = interp1([a1(n) a2(n)],[f1(n) f2(n)],x(ind));
% 		else % Injection is within sampling points
% 			ind1 = find(x<=a2(n));
% 			ind2 = find(x>=a1(n));
% 			out(ind1) = interp1([a1(n) a2(n)+C],[f1(n) f2(n)],C+x(ind1));
% 			out(ind2) = interp1([a1(n) a2(n)+C],[f1(n) f2(n)],x(ind2));
% 		end
% 	end
% end

