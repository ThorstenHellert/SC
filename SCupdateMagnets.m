function SC = SCupdateMagnets(SC,varargin)
% SCupdateMagnets
% ===============
%
% NAME
% ----
% SCupdateMagnets - Updates the magnetic fields in RING
%
% SYNOPSIS
% --------
% `SC = SCupdateMagnets(SC [, ords])`
%
%
% DESCRIPTION
% -----------
% Updates the magnets specified in `SC.RING` as specified in `ords`. If no ordinates are given
% explicitly, all registered magnets defined in `SC.ORD.Magnet` are updated. For each magnet the
% setpoints (`SetPointA/B`) and calibration errors (`CalErrorA/B`) are evaluated.
% If systematic multipole components are specified, e.g. in `SysPolBFromB` for systematic 
% PolynomB-multipoles induced by PolynomB entries, the corresponding multipole components are scaled 
% by the current magnet excitation and added, as well as static field offsets (if specified in 
% `PolynomA/BOffset`). 
% If the considered magnet has a bending angle error (from pure bending angle error or due to a
% combined function magnet), the corresponding horizontal dipole magnetic field is calculated and
% added to the PolynomB(1) term. It is thereby assured that a dipole error doesn't alter the
% coordinate system.
%
% If the considered magnet is registered as a split magnet (`'MasterOf'`), the errors and setpoints
% of the master magnet are applied to the fields of the child magnets. Note that split quadrupole
% magnets with different gradients, however, or split CMs can currently not be updated correctly.
%
% RETURN VALUE
% ------------
% `SC`::
% 	Base structure with updated magnets.
%
% SEE ALSO
% --------
% *SCregisterMagnets*, *SCapplyErrors*, *SCsetMultipoles*, *SCsetMags2SetPoints*, *SCsetCMs2SetPoints*


	% Check if ordinates are given explicitly. If not, take all registered magnets
	if isempty(varargin)
		ords = SC.ORD.Magnet;
	else
		ords = varargin{1};
	end

	% Loop over all specified magnets
	for ord=ords
		SC = updateMagnets(SC,ord,ord);

		% Check for Master/Child magnets and update children
		if isfield(SC.RING{ord},'MasterOf')
			for childOrd=SC.RING{ord}.MasterOf
				SC = updateMagnets(SC,ord,childOrd);
			end
		end
	end
end


% Update the magnetic fields
function SC = updateMagnets(SC,source,target)
	% Writes the magnetic fields in magnet specified in 'target' by reading the magnet setting and
	% errors of the magnet defined in 'source'.

	% Apply magnet setpoint including calibration error
	SC.RING{target}.PolynomB = SC.RING{source}.SetPointB .* addPadded(ones(size(SC.RING{source}.SetPointB)),SC.RING{source}.CalErrorB);
	SC.RING{target}.PolynomA = SC.RING{source}.SetPointA .* addPadded(ones(size(SC.RING{source}.SetPointA)),SC.RING{source}.CalErrorA);

	% Initialize temporary systematic multipole arrays
	sysPolynomB = [];
	sysPolynomA = [];
	
	% Check for any systematic PolynomB multipole errors from PolynomB entrys
	if isfield(SC.RING{target},'SysPolBFromB')
		for n=1:length(SC.RING{target}.SysPolBFromB)
			if ~isempty(SC.RING{target}.SysPolBFromB{n})
				% Scale PolynomB multipoles by current PolynomB magnet excitation
				sysPolynomB{end+1} = SC.RING{target}.PolynomB(n) * SC.RING{target}.SysPolBFromB{n};
			end
		end
	end
	% Check for any systematic PolynomB multipole errors from PolynomA entrys
	if isfield(SC.RING{target},'SysPolBFromA')
		for n=1:length(SC.RING{target}.SysPolBFromA)
			if ~isempty(SC.RING{target}.SysPolBFromA{n})
				% Scale PolynomB multipoles by current PolynomA magnet excitation
				sysPolynomB{end+1} = SC.RING{target}.PolynomA(n) * SC.RING{target}.SysPolBFromA{n};
			end
		end
	end
	
	% Check for any systematic PolynomA multipole errors from PolynomB entrys
	if isfield(SC.RING{target},'SysPolAFromB')
		for n=1:length(SC.RING{target}.SysPolAFromB)
			if ~isempty(SC.RING{target}.SysPolAFromB{n})
				% Scale PolynomA multipoles by current PolynomB magnet excitation
				sysPolynomA{end+1} = SC.RING{target}.PolynomB(n) * SC.RING{target}.SysPolAFromB{n};
			end
		end
	end
	% Check for any systematic PolynomA multipole errors from PolynomA entrys
	if isfield(SC.RING{target},'SysPolAFromA')
		for n=1:length(SC.RING{target}.SysPolAFromA)
			if ~isempty(SC.RING{target}.SysPolAFromA{n})
				% Scale PolynomA multipoles by current PolynomA magnet excitation
				sysPolynomA{end+1} = SC.RING{target}.PolynomA(n) * SC.RING{target}.SysPolAFromA{n};
			end
		end
	end
	
	
	% Add systematic multipoles for PolynomA
	if ~isempty(sysPolynomA)
		% Sum up all individual contributions...
		for n=1:length(sysPolynomA)-1
			sysPolynomA{n+1} = addPadded(sysPolynomA{n+1},sysPolynomA{n});
		end
		% ...and add them to PolynomA field
		SC.RING{target}.PolynomA = addPadded(SC.RING{target}.PolynomA, sysPolynomA{end});
	end
	% Add systematic multipoles for PolynomB
	if ~isempty(sysPolynomB)
		% Sum up all individual contributions...
		for n=1:length(sysPolynomB)-1
			sysPolynomB{n+1} = addPadded(sysPolynomB{n+1},sysPolynomB{n});
		end
		% ...and add them to PolynomB field
		SC.RING{target}.PolynomB = addPadded(SC.RING{target}.PolynomB, sysPolynomB{end});
	end
	
		
	% Add static multipole errors of target magnet
	SC.RING{target}.PolynomB = addPadded(SC.RING{target}.PolynomB,SC.RING{target}.PolynomBOffset);
	SC.RING{target}.PolynomA = addPadded(SC.RING{target}.PolynomA,SC.RING{target}.PolynomAOffset);
	

	% Add bending angle errors
	if isfield(SC.RING{source},'BendingAngleError')
		SC.RING{target}.PolynomB(1) = SC.RING{target}.PolynomB(1) + SC.RING{source}.BendingAngleError * SC.RING{target}.BendingAngle/SC.RING{target}.Length;
	end

	% Include bending angle error from combinded function magnets
	if isfield(SC.RING{source},'BendingAngle')
		% Check if magnet is registered as combined function (bending angle depends on quad component)
		if isfield(SC.RING{source},'CombinedFunction') && SC.RING{source}.CombinedFunction==1
			% Compute quadrupole deviation from design value
			alpha_act = SC.RING{source}.SetPointB(2) * (1+SC.RING{source}.CalErrorB(2)) / SC.RING{source}.NomPolynomB(2);
			% Compute effective bending angle
			effBendingAngle = alpha_act * SC.RING{target}.BendingAngle;
		
			% Add bending angle difference to PolynomB field
			SC.RING{target}.PolynomB(1) = SC.RING{target}.PolynomB(1) + (effBendingAngle - SC.RING{target}.BendingAngle) / SC.RING{target}.Length;
		end
	end

	

	% Deal with CorrectorPass
	if strcmp(SC.RING{source}.PassMethod,'CorrectorPass')
		SC.RING{target}.KickAngle(1) = SC.RING{target}.PolynomB(1);
		SC.RING{target}.KickAngle(2) = SC.RING{target}.PolynomA(1);
	end
	
	% Adjust order for tracking
	SC.RING{target}.MaxOrder=length(SC.RING{target}.PolynomB)-1;

end



% Auziliary function to add unequally long arrays (zero padding)
function v = addPadded(v1,v2)
	if ~((iscolumn(v1)&&iscolumn(v2))||(isrow(v1)&&isrow(v2)))
		error('Wrong dimensions.');
	end
	l1=length(v1);
	l2=length(v2);
	if l2>l1; v1(l2)=0; end
	if l1>l2; v2(l1)=0; end
	v=v1+v2;
end
