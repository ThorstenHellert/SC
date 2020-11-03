function setpoints = SCgetCMSetPoints(SC,CMords,nDim)
% SCgetCMSetPoints
% ================
%
% NAME
% ----
% SCgetCMSetPoints - Return current CM setpoints
%
% SYNOPSIS
% --------
% `setpoints = SCgetCMSetPoints(SC, CMords, nDim)`
%
%
% DESCRIPTION
% -----------
% Reads the setpoints of the CMs specified in `CMords` in the dimension `nDim`.
%
% INPUTS
% ------
% `SC`::
%	SC base structure
% `CMords`::
%	Array of CM ordinates in the lattice structure
% `nDim`::
%	Integer specifying CM dimension ([1|2] -> [hor|ver])
%
%
% RETURN VALUES
% -------------
% `setpoints`::
%	CM setpoints [rad]
%
% SEE ALSO
% --------
% *SCsetCMs2SetPoints*, *SCregisterMagnets*


	setpoints = nan(1,length(CMords));

	% Set setpoints and update magnets
	for idx = 1:length(CMords)

		% Get setpoint normalization factor (field <-> kick)
		if strcmp(SC.RING{CMords(idx)}.PassMethod,'CorrectorPass')
			normBy = [1 1];
		else
			normBy = [-1 1]*SC.RING{CMords(idx)}.Length; % positive setpoint -> positive kick -> negative horizontal field
		end

		% Read setpoints
		if nDim==1
			setpoints(idx) = SC.RING{CMords(idx)}.SetPointB(1) * normBy(nDim);
		else
			setpoints(idx) = SC.RING{CMords(idx)}.SetPointA(1) * normBy(nDim);
		end

	end
end

