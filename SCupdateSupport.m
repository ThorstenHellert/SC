function SC = SCupdateSupport(SC,varargin)
% SCupdateSupport
% ===============
%
% NAME
% ----
% SCupdateSupport - updates the misalignments resulting from the support structure
%
% SYNOPSIS
% --------
% `SC = SCupdateSupport(SC [, options])`
%
%
% DESCRIPTION
% -----------
% This function updates the offsets and rolls of the elements in `SC.RING`
% based on the current support errors, by setting the lattice fields `T1`, `T2`, and
% `R1`, `R2` for magnets and the fields `SupportOffset` and `SupportRoll` for BPMs.
%
%
% INPUT
% -----
% `SC`::
%	Initial `SC` base structure.
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'BPMstructOffset'` (`1`)::
%	If true, BPM offsets are updated.
% `'MAGstructOffset'` (`1`)::
%	If true, magnet offsets are updated.
%
%
% RETURN VALUE
% ------------
% `SC` ::
%	Base structure with updated `SC.RING`.
%
% SEE ALSO
% --------
% *SCregisterSupport*, *SCgetSupportOffset*, *SCgetSupportRoll*, *SCplotSupport*


	p = inputParser;
	p.KeepUnmatched=true;
	addParameter(p,'BPMstructOffset',1);
	addParameter(p,'MAGstructOffset',1);
	parse(p,varargin{:});

	if p.Results.MAGstructOffset
		% Process Magnets
		if isfield(SC.ORD,'Magnet')
			ords    = unique([SC.ORD.Magnet]);
			s       = findspos(SC.RING,ords);
			offsets = SCgetSupportOffset(SC,s);
			rolls   = SCgetSupportRoll(SC,s);
			
			for i=1:length(ords)
				ord = ords(i);

				% Update support offset
				SC.RING{ord}.SupportOffset = offsets(:,i)';

				% Update support roll angles
				SC.RING{ord}.SupportRoll = rolls(:,i)';

				% Get length of magnet 
				magLength   = SC.RING{ord}.Length;
				
				% Get bending angle of magnet 
				if isfield(SC.RING{ord},'BendingAngle')
					magTheta = SC.RING{ord}.BendingAngle;
				else
					magTheta = 0;
				end
				
				% Magnet horizontal, vertical and longitudinal offset
				dx = SC.RING{ord}.SupportOffset(1) + SC.RING{ord}.MagnetOffset(1);
				dy = SC.RING{ord}.SupportOffset(2) + SC.RING{ord}.MagnetOffset(2);
				dz = SC.RING{ord}.SupportOffset(3) + SC.RING{ord}.MagnetOffset(3);
				% Magnet roll around z-, x- and yaxis
				az = SC.RING{ord}.MagnetRoll(1) + SC.RING{ord}.SupportRoll(1);
				ax = SC.RING{ord}.MagnetRoll(2) + SC.RING{ord}.SupportRoll(2);
				ay = SC.RING{ord}.MagnetRoll(3) + SC.RING{ord}.SupportRoll(3);

				% Calculate magnet transformations 
				[T1,T2,R1,R2] = SCgetTransformation(dx,dy,dz,ax,ay,az,magTheta,magLength);
				
				% Set fields in lattice element
				SC.RING{ord}.T1 = T1;
				SC.RING{ord}.T2 = T2;
				SC.RING{ord}.R1 = R1;
				SC.RING{ord}.R2 = R2;
				
				
				% Check for Master/Child magnets and copy error
				if isfield(SC.RING{ord},'MasterOf')
					for childOrd=SC.RING{ord}.MasterOf
						SC.RING{childOrd}.T1=SC.RING{ord}.T1;
						SC.RING{childOrd}.T2=SC.RING{ord}.T2;
						% TODO: Adjust following lines to properly manage split longitudinal gradient dipoles
						SC.RING{childOrd}.R1=SC.RING{ord}.R1;
						SC.RING{childOrd}.R2=SC.RING{ord}.R2;
					end
				end
			end
		else
			warning('SC: No magnets have been registered!')
		end
	end

	
	if p.Results.BPMstructOffset
		% Process BPMs`
		if isfield(SC.ORD,'BPM')
			ords    = unique([SC.ORD.BPM]);
			s       = findspos(SC.RING,ords);
			offsets = SCgetSupportOffset(SC,s);
			rolls   = SCgetSupportRoll(SC,s);
			
			
			for i=1:length(ords)
				ord = ords(i);
			
				% Update support offset
				SC.RING{ord}.SupportOffset = offsets(1:2,i)'; % Longitudinal BPM offsets not yet implemented

				% Update support roll angles
				SC.RING{ord}.SupportRoll = rolls(1,i)'; % BPM pitch and yaw angles not yet implemented
		
			end
		else
			warning('SC: No BPMs have been registered!')
		end
	end		
					
end
