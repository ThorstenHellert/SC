function SC = SCapplyErrors(SC,varargin)
% SCapplyErrors
% =============
%
% NAME
% ----
% SCapplyErrors - Apply errors to lattice and diagnostic devices
%
% SYNOPSIS
% --------
% `SC = SCapplyErrors(SC)`
%
% INPUT
% -----
% `SC`::
% 	The SC base structure.
%
% DESCRIPTION
% -----------
% Applies errors to cavities, injection trajectory, BPMs, circumference,
% support structures and magnets if the corresponding uncertanties defined in
% `SC.SIG` are set. For example, for a magnet with ordinate `ord` every field
% defined in `SC.SIG.Mag{ord}` will be used to generate a random number using a
% Gaussian distribution with a cutoff at `+/-2sigma` and `sigma` being the
% value of the uncertainty field. The number will be stored in the
% corresponding field of the lattice structure, thus `SC.RING{ord}`. An
% exeption are bending angle errors which are stored in the `BendingAngleError`
% field. See examples in the register functions for more details.
%
% RETURN VALUE
% ------------
% `SC`::
% 	The SC base structure with applied errors.
%
% SEE ALSO
% --------
% *SCregisterMagnets*, *SCregisterSupport*, *SCregisterBPMs*, *SCregisterCAVs*, *SCrampUpErrors*

	if ~isfield(SC,'SIG')
		return;
	end

	% Apply cavity errors
	SC = applyCavityError(SC);

	% Apply injected beam errors
	SC = applyInjectionError(SC);

	% Apply BPM errors
	SC = applyBPMerrors(SC);

	% Apply circumference error
	SC = applyCircumferenceError(SC);

	% Apply magnet support errors
	SC = applySupportAlignmentError(SC);

	% Apply magnet errors
	SC = applyMagnetError(SC);

	% Update magnet support model
	SC = SCupdateSupport(SC,varargin{:});

	% Update magnetic fields
	if isfield(SC.ORD,'Magnet')
		SC = SCupdateMagnets(SC);
	end
	
	% Update cavity fields
	if isfield(SC.ORD,'Cavity')
		SC = SCupdateCAVs(SC);
	end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cavity errors
function SC = applyCavityError(SC)
	if isfield(SC.SIG,'RF')
		for ord=SC.ORD.Cavity
			
			% Check if uncertainty for current BPM is given
			if isempty(SC.SIG.RF{ord})
				continue
			end

			% Loop over uncertainties
			for field=fieldnames(SC.SIG.RF{ord})'
				SC.RING{ord}.(field{1}) = SC.SIG.RF{ord}.(field{1}) * SCrandnc(2);
			end
		end
	end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Injected beam errors
function SC = applyInjectionError(SC)
	if isfield(SC.SIG,'staticInjectionZ')
		% Apply systematic injection error
		SC.INJ.Z0               = SC.INJ.Z0ideal + SC.SIG.staticInjectionZ(:) .* SCrandnc(2,6,1);
		fprintf('Static injection error applied.\n');
	end
	if isfield(SC.SIG,'randomInjectionZ')
		% Define random injection error (shot-to-shot)
		SC.INJ.randomInjectionZ = SC.SIG.randomInjectionZ(:);
		fprintf('Random injection error applied.\n');
	end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BPM errors
function SC = applyBPMerrors(SC)
	% Check if any uncertainties for the magnets are given
	if ~isfield(SC.SIG,'BPM')
		return;
	end
	
	% Loop over BPMs
	for ord = SC.ORD.BPM
	
		% Check if uncertainty for current BPM is given
		if isempty(SC.SIG.BPM{ord})
			continue
		end
	
		% Loop over uncertainties
		for field=fieldnames(SC.SIG.BPM{ord})'
			
			% Noise value should be copied into final structure, all others randomly generated
			if regexp(field{1},'Noise')
				SC.RING{ord}.(field{1}) = SC.SIG.BPM{ord}.(field{1});
			else
				SC.RING{ord}.(field{1}) = SC.SIG.BPM{ord}.(field{1}) .* SCrandnc(2,size(SC.SIG.BPM{ord}.(field{1})));
			end
		end
	end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Circumference error
function SC = applyCircumferenceError(SC)
	if isfield(SC.SIG,'Circumference')
		% Define circumference error
		circScaling = 1 + SC.SIG.Circumference * SCrandnc(2,1,1);
		% Apply circumference error
		SC.RING = SCscaleCircumference(SC.RING,circScaling,'rel');
		fprintf('Circumference error applied.\n');
	end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suport structure alignemt errors
function SC = applySupportAlignmentError(SC)

	% Loop over different support types
	for type = {'Girder','Plinth','Section'}
		% Check if support type is registered
		if ~isfield(SC.ORD,type{1})
			continue;
		end

		% Loop over support structure pairs
		for ordPair=SC.ORD.(type{1})
			% Check if uncertainties are set
			if isempty(SC.SIG.Support{ordPair(1)})
				continue;
			end

			% Loop over fields in sigma structure
			for field=fieldnames(SC.SIG.Support{ordPair(1)})'
				% Check if considered field in sigma structure is related to currently evaluated support structure type
				if ~regexp(field{1},type{1})
					continue;
				end

				% Generate random error for support structure beginning
				SC.RING{ordPair(1)}.(field{1}) = SC.SIG.Support{ordPair(1)}.(field{1}) .* SCrandnc(2,size(SC.SIG.Support{ordPair(1)}.(field{1})));

				% Check if uncertanty is specified for endpoint
				if length(SC.SIG.Support)>=ordPair(2) && isfield(SC.SIG.Support{ordPair(2)},field{1})
					% Generate random error for support structure endpoint
					SC.RING{ordPair(2)}.(field{1}) = SC.SIG.Support{ordPair(2)}.(field{1}) .* SCrandnc(2,size(SC.SIG.Support{ordPair(2)}.(field{1})));
				else
					% Copy support structure endpoint from structure beginning
					SC.RING{ordPair(2)}.(field{1}) = SC.RING{ordPair(1)}.(field{1});
				end
			end
		end
	end
	
	% Get girder pitch and yaw angles
	if isfield(SC.ORD,'Girder')
		% Loop over girders
		for ordPair=SC.ORD.Girder
			% Get girder length
			gLength = abs(diff(findspos(SC.RING,ordPair)));
			
			% Check if pitch and yaw errors are given explicitly
			if length(SC.RING{ordPair(1)}.GirderRoll)==3 && (SC.RING{ordPair(1)}.GirderRoll(2)~=0 || SC.RING{ordPair(1)}.GirderRoll(3)~=0)
				% Check if girder start- and endpoints have individual offset uncertainties
				if length(SC.SIG.Support)>=ordPair(2) && ~isempty(SC.SIG.Support{ordPair(2)}) && isfield(SC.SIG.Support{ordPair(2)},'GirderOffset')
					error('Pitch or yaw angle errors can not be given explicitly if girder start and endpoints each have offset uncertainties.')
				end
				
				% Adjust girder startpoint vertcial offset
				SC.RING{ordPair(1)}.GirderOffset(2) = SC.RING{ordPair(1)}.GirderOffset(2) - SC.RING{ordPair(1)}.GirderRoll(2)*gLength;
				% Adjust girder endpoint vertcial offset
				SC.RING{ordPair(2)}.GirderOffset(2) = SC.RING{ordPair(2)}.GirderOffset(2) + SC.RING{ordPair(1)}.GirderRoll(2)*gLength;
				% Adjust girder startpoint horizontal offset
				SC.RING{ordPair(1)}.GirderOffset(1) = SC.RING{ordPair(1)}.GirderOffset(1) - SC.RING{ordPair(1)}.GirderRoll(3)*gLength;
				% Adjust girder endpoint horizontal offset
				SC.RING{ordPair(2)}.GirderOffset(1) = SC.RING{ordPair(2)}.GirderOffset(1) + SC.RING{ordPair(1)}.GirderRoll(3)*gLength;
				
			else
				% Get girder pitch from horizontal start- and endpoint offsets
				SC.RING{ordPair(1)}.GirderRoll(2) = (SC.RING{ordPair(2)}.GirderOffset(1) - SC.RING{ordPair(1)}.GirderOffset(1))/gLength;
				% Get girder yaw from vertical start- and endpoint offsets
				SC.RING{ordPair(1)}.GirderRoll(3) = (SC.RING{ordPair(2)}.GirderOffset(2) - SC.RING{ordPair(1)}.GirderOffset(2))/gLength;
				
			end
		end
	end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnet errors
function SC = applyMagnetError(SC)
	% Check if any uncertainties for the magnets are given
	if ~isfield(SC.SIG,'Mag')
		return;
	end
	for ord = SC.ORD.Magnet

		if isempty(SC.SIG.Mag{ord})
			continue
		end
		for field=fieldnames(SC.SIG.Mag{ord})'
			% Bending angle error gets applied differently
			if strcmp(field{1},'BendingAngle')
				SC.RING{ord}.BendingAngleError = SC.SIG.Mag{ord}.BendingAngle * SCrandnc(2,1,1);
			else
				SC.RING{ord}.(field{1}) = SC.SIG.Mag{ord}.(field{1}) .* SCrandnc(2,size(SC.SIG.Mag{ord}.(field{1})));
			end
		end
	end
end


