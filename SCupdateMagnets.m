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
% Updates the magnets specified in `RING` as specified in `ords`. If no ordinates are given
% explicitly, all registered magnets defined in `SC.ORD.Magnet` are updated. For each magnet the
% setpoints (`SetPointA/B`) and calibration errors (`CalErrorA/B`) are evaluated and field offsets
% (`PolynomA/BOffset`) are added.
% If the considered magnet has a bending angle, the corresponding horizontal dipole magnetic field
% is calculated. Thereafter, the fields induced by a roll error are calculated. Finally the fields
% described by the design bending angle (if present) are subtracted from the `PolynomA/B` terms. It is
% thereby assured that a rotation of dipole magnet doesn't alter the coordinate system.
%
% If the considered magnet is registered as a slpit magnet (`'MasterOf'`), the errors and setpoints 
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

		% Check for Master/Child magnets and update childs
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

	% Add bending angle errors
	if isfield(SC.RING{source},'BendingAngleError')
		SC.RING{target}.PolynomB(1) = SC.RING{target}.PolynomB(1) + SC.RING{source}.BendingAngleError * SC.RING{target}.BendingAngle/SC.RING{target}.Length;
	end

	% Check for multipole errors of 'write' magnet
	if isfield(SC.RING{target},'PolynomBOffset')
		SC.RING{target}.PolynomB = addPadded(SC.RING{target}.PolynomB,SC.RING{target}.PolynomBOffset);
		SC.RING{target}.PolynomA = addPadded(SC.RING{target}.PolynomA,SC.RING{target}.PolynomAOffset);
	end

	% Include bending angle temporarily in PolynomB field
	if isfield(SC.RING{source},'BendingAngle')
		% Check if magnet is registered as combined function (bending angle depends on quad component)
		if isfield(SC.RING{source},'CombinedFunction') && SC.RING{source}.CombinedFunction==1
			% Compute quadrupole deviation from design value
			alpha_act = SC.RING{source}.SetPointB(2) * (1+SC.RING{source}.CalErrorB(2)) / SC.RING{source}.NomPolynomB(2);
			% Compute effective bending angle
			effBendingAngle = alpha_act * SC.RING{target}.BendingAngle;
		else
			effBendingAngle = SC.RING{target}.BendingAngle;
		end
		% Add bending angle to field
		SC.RING{target}.PolynomB(1) = SC.RING{target}.PolynomB(1) + effBendingAngle / SC.RING{target}.Length;
	end

	% Apply roll errors. See Accelerator Handbook p445.
	C = cos(SC.RING{source}.RollAngle * [1:length(SC.RING{target}.PolynomB)]);
	S = sin(SC.RING{source}.RollAngle * [1:length(SC.RING{target}.PolynomB)]);
	PolynomBprime = SC.RING{target}.PolynomB .* C - SC.RING{target}.PolynomA .* S;
	PolynomAprime = SC.RING{target}.PolynomB .* S + SC.RING{target}.PolynomA .* C;
	SC.RING{target}.PolynomB = PolynomBprime;
	SC.RING{target}.PolynomA = PolynomAprime;

	% Remove bending angle from PolynomB field
	if isfield(SC.RING{target},'BendingAngle')
		SC.RING{target}.PolynomB(1) = SC.RING{target}.PolynomB(1) - SC.RING{target}.BendingAngle / SC.RING{target}.Length;
	end

	% Deal with CorrectorPass
	if strcmp(SC.RING{source}.PassMethod,'CorrectorPass')
		SC.RING{target}.KickAngle(1) = SC.RING{target}.PolynomB(1);
		SC.RING{target}.KickAngle(2) = SC.RING{target}.PolynomA(1);
	end
	
	% Increase order for tracking
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
