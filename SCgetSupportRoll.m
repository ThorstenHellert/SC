function roll = SCgetSupportRoll(SC,s)
% SCgetSupportRoll
% ================
%
% NAME
% ----
% SCgetSupportRoll - Calculates the combined support structure roll, pitch and yaw angles at a certain location
%
% SYNOPSIS
% --------
% `roll = SCgetSupportRoll(SC, s)`
%
%
% DESCRIPTION
% -----------
% This function evaluates the total roll, pitch and yaw angles of the support structures that have
% been defined via *SCregisterSupport* at the longitudinal positions `s`.
% The support structure pitch and yaw angles are calculated from the horizontal and vertical offsets
% of the start end endpoints of the corresponding top layer support structure 
% (the order is: girders->plinths->sections). The roll angle is a sum of the roll angles of all
% underlying support structures.
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
% `roll`::
%	`[3,length(s)]`-array containing the [az/ax/ay] total support structure rolls at `s`.
%
% SEE ALSO
% --------
% *SCregisterSupport*, *SCupdateSupport*, *SCgetSupportOffset*, *SCplotSupport*, *SCapplyErrors*

	% s-positions along the ring
	s0 = findspos(SC.RING,1:length(SC.RING));

	% Initialize rolls along the ring
	roll0 = zeros(3,length(s0));
	
	% Read elements length from RING
	for n=1:length(SC.RING)
		lengths(n) = SC.RING{n}.Length;
	end
	
	% Circumference
	C = sum(lengths);    

	% Read support structure offset along the ring
	off0 = SCgetSupportOffset(SC,s0);
	
	% Get start and end ordinates of all support structure elements
	supportOrds = getSupportOrds(SC);
	
	% Loop over support structure elements ('Section'->'Plinth'->'Girder')
	for type = {'Section','Plinth','Girder'}
		% Check if support structure is registered
		if isfield(supportOrds,type{1})
			for nEl=1:length(supportOrds.(type{1}))
				ords = supportOrds.(type{1}){nEl};
				if diff(ords)>0
					% Add roll error from current support structure
					roll0(1,ords(1):ords(2)) = roll0(1,ords(1):ords(2)) + SC.RING{ords(1)}.([type{1} 'Roll'])(1);
					% Overwrite pitch and yaw angles from current support structure hor. and ver. offsets
					roll0(2,ords(1):ords(2)) = (off0(2,ords(2))-off0(2,ords(1))) / (s0(ords(2))-s0(ords(1)));
					roll0(3,ords(1):ords(2)) = (off0(1,ords(2))-off0(1,ords(1))) / (s0(ords(2))-s0(ords(1)));
				else
					% Add roll error from current support structure
					roll0(1,ords(1):end) = roll0(1,ords(1):end) + SC.RING{ords(1)}.([type{1} 'Roll'])(1);
					roll0(1,1:ords(2))   = roll0(1,1:ords(2))   + SC.RING{ords(1)}.([type{1} 'Roll'])(1);
					% Overwrite pitch angle from current support structure vertical offsets
					roll0(2,ords(1):end) = (off0(2,ords(2))-off0(2,ords(1))) / (C-s0(ords(1))-s0(ords(2)));
					roll0(2,1:ords(2))   = (off0(2,ords(2))-off0(2,ords(1))) / (C-s0(ords(1))-s0(ords(2)));
					% Overwrite yaw angle from current support structure horiaontal offsets
					roll0(3,ords(1):end) = (off0(1,ords(2))-off0(1,ords(1))) / (C-s0(ords(1))-s0(ords(2)));
					roll0(3,1:ords(2))   = (off0(1,ords(2))-off0(1,ords(1))) / (C-s0(ords(1))-s0(ords(2)));
				end
			end
		end
	end

	% Get rolls at requested s-positions
	[~,b] = unique(s0);
	roll(1,:) = interp1(s0(b),roll0(1,b),s,'linear','extrap');
	roll(2,:) = interp1(s0(b),roll0(2,b),s,'linear','extrap');
	roll(3,:) = interp1(s0(b),roll0(3,b),s,'linear','extrap');

end

function supportOrds = getSupportOrds(SC)
	
	supportOrds = [];

	% Loop over support structure types
	for type = {'Section','Plinth','Girder'}
		% Check if support structure is registered
		if isfield(SC.ORD,type{1})
			for i=1:size(SC.ORD.(type{1}),2)
				supportOrds.(type{1}){i} = SC.ORD.(type{1})(:,i);
			end
		end
	end
end

