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
% `RollAngle` in the elements.
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

% TODO: Centralize `MasterOf` functionality

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

			for i=1:length(ords)
				ord = ords(i);
				off = offsets(:,i)';

				% Update support offset
				SC.RING{ord}.SupportOffset = off;

				% Update magnet offset
				T = zeros(6,1);
				T([1,3]) = SC.RING{ord}.SupportOffset + SC.RING{ord}.MagnetOffset;
				SC.RING{ord}.T1 = -T;
				SC.RING{ord}.T2 = +T;

				% Check for Master/Child magnets and copy error
				if isfield(SC.RING{ord},'MasterOf')
					for childOrd=SC.RING{ord}.MasterOf
						SC.RING{childOrd}.T1=SC.RING{ord}.T1;
						SC.RING{childOrd}.T2=SC.RING{ord}.T2;
					end
				end
				
				% Get girder rolls for magnets
				if isfield(SC.ORD,'Girder')
					% Find girder index of magnet
					gInd = intersect(find(ord>SC.ORD.Girder(1,:)),find(ord<SC.ORD.Girder(2,:)));
					if ~isempty(gInd)
						% Write girder roll in magnet element
						SC.RING{ord}.SupportRoll = SC.RING{SC.ORD.Girder(1,gInd)}.GirderRoll;
					end
				end
				
				% Update overall magnet roll error
				SC.RING{ord}.RollAngle = SC.RING{ord}.MagnetRoll + SC.RING{ord}.SupportRoll;
				
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
			
			
			for i=1:length(ords)
				ord = ords(i);
				off = offsets(:,i)';
			
				% Update support offset
				SC.RING{ord}.SupportOffset = off;

				% Get girder rolls for BPMs
				if isfield(SC.ORD,'Girder')
					% Find girder index of BPM
					gInd = intersect(find(ord>SC.ORD.Girder(1,:)),find(ord<SC.ORD.Girder(2,:)));
					if ~isempty(gInd)
						% Write girder roll in BPM element
						SC.RING{ord}.SupportRoll = SC.RING{SC.ORD.Girder(1,gInd)}.GirderRoll;
					end
				end
				
			end
		else
			warning('SC: No BPMs have been registered!')
		end
	end		
					
end
