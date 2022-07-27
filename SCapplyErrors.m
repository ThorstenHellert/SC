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
% `SC = SCapplyErrors(SC [, options])`
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
% Gaussian distribution with a cutoff (see option below) and `sigma` being the
% value of the uncertainty field. The number will be stored in the
% corresponding field of the lattice structure, thus `SC.RING{ord}`. An
% exeption are bending angle errors which are stored in the `BendingAngleError`
% field. See examples in the SCregister* functions for more details.
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'nSig'` (2)::
%	Number of sigmas at which the Gaussian distribution of errors is truncated if not defined
%   explicitly for individual uncertainties.
%
% RETURN VALUE
% ------------
% `SC`::
% 	The SC base structure with applied errors.
%
% SEE ALSO
% --------
% *SCregisterMagnets*, *SCregisterSupport*, *SCregisterBPMs*, *SCregisterCAVs*, *SCrampUpErrors*

	% Parse optional arguments
	p = inputParser;
	addOptional(p,'nSig',2);
	parse(p,varargin{:});
	par=p.Results;

	if ~isfield(SC,'SIG')
		warning('No uncertanties provided.')
		return;
	end

	
	
	% Apply cavity errors
	SC = applyCavityError(SC,par);

	% Apply injected beam errors
	SC = applyInjectionError(SC,par);

	% Apply BPM errors
	SC = applyBPMerrors(SC,par);

	% Apply circumference error
	SC = applyCircumferenceError(SC,par);

	% Apply magnet support errors
	SC = applySupportAlignmentError(SC,par);

	% Apply magnet errors
	SC = applyMagnetError(SC,par);

	% Update magnet support model
	SC = SCupdateSupport(SC);

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
function SC = applyCavityError(SC,par)
	if isfield(SC.SIG,'RF')
		for ord=SC.ORD.Cavity
			
			% Check if uncertainty for current cavity is given
			if isempty(SC.SIG.RF{ord})
				continue
			end

			% Loop over uncertainties
			for field=fieldnames(SC.SIG.RF{ord})'
				% Check if individual cutoff is defined
				if iscell(SC.SIG.RF{ord}.(field{1}))
					SC.RING{ord}.(field{1}) = SC.SIG.RF{ord}.(field{1}){1} * SCrandnc(SC.SIG.RF{ord}.(field{1}){2});
				else
					SC.RING{ord}.(field{1}) = SC.SIG.RF{ord}.(field{1}) * SCrandnc(par.nSig);
				end
			end
		end
	end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Injected beam errors
function SC = applyInjectionError(SC,par)
	if isfield(SC.SIG,'staticInjectionZ')
		% Apply systematic injection error
		SC.INJ.Z0               = SC.INJ.Z0ideal + SC.SIG.staticInjectionZ(:) .* SCrandnc(par.nSig,6,1);
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
function SC = applyBPMerrors(SC,par)
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
				% Check if individual cutoff is defined
				if iscell(SC.SIG.BPM{ord}.(field{1}))
					SC.RING{ord}.(field{1}) = SC.SIG.BPM{ord}.(field{1}){1} .* SCrandnc(SC.SIG.BPM{ord}.(field{1}){2},size(SC.SIG.BPM{ord}.(field{1}){1}));
				else
					SC.RING{ord}.(field{1}) = SC.SIG.BPM{ord}.(field{1}) .* SCrandnc(par.nSig,size(SC.SIG.BPM{ord}.(field{1})));
				end
			end
		end
	end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Circumference error
function SC = applyCircumferenceError(SC,par)
	if isfield(SC.SIG,'Circumference')
		% Define circumference error
		circScaling = 1 + SC.SIG.Circumference * SCrandnc(par.nSig,1,1);
		% Apply circumference error
		SC.RING = SCscaleCircumference(SC.RING,circScaling,'rel');
		fprintf('Circumference error applied.\n');
	end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suport structure alignemt errors
function SC = applySupportAlignmentError(SC,par)
	
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
				if isempty(strfind(field{1},type{1}))
					continue;
				end

				% Generate random error for support structure beginning
				if iscell(SC.SIG.Support{ordPair(1)}.(field{1}))
					% Individual cutoff is given
					SC.RING{ordPair(1)}.(field{1}) = SC.SIG.Support{ordPair(1)}.(field{1}){1} .* SCrandnc(SC.SIG.Support{ordPair(1)}.(field{1}){2},size(SC.SIG.Support{ordPair(1)}.(field{1}){1}));
				else
					SC.RING{ordPair(1)}.(field{1}) = SC.SIG.Support{ordPair(1)}.(field{1}) .* SCrandnc(par.nSig,size(SC.SIG.Support{ordPair(1)}.(field{1})));
				end

				% Check if uncertanty is specified for endpoint
				if length(SC.SIG.Support)>=ordPair(2) && isfield(SC.SIG.Support{ordPair(2)},field{1})
					% Generate random error for support structure endpoint
					if iscell(SC.SIG.Support{ordPair(2)}.(field{1}))
						% Individual cutoff is given
						SC.RING{ordPair(2)}.(field{1}) = SC.SIG.Support{ordPair(2)}.(field{1}){1} .* SCrandnc(SC.SIG.Support{ordPair(2)}.(field{1}){2},size(SC.SIG.Support{ordPair(2)}.(field{1}){1}));
					else
						SC.RING{ordPair(2)}.(field{1}) = SC.SIG.Support{ordPair(2)}.(field{1}) .* SCrandnc(par.nSig,size(SC.SIG.Support{ordPair(2)}.(field{1})));
					end
				else
					% Copy support structure endpoint from structure beginning
					SC.RING{ordPair(2)}.(field{1}) = SC.RING{ordPair(1)}.(field{1});
				end
			end

			% Get support structure length
			if diff(ordPair)>=0
				structLength = abs(diff(findspos(SC.RING,ordPair)));
			else
				structLength = findspos(SC.RING,ordPair(2)) + diff(findspos(SC.RING,[ordPair(1) length(SC.RING)+1]));
			end
			
			% Check if support structure pitch errors are given explicitly
			if SC.RING{ordPair(1)}.([type{1} 'Roll'])(2)~=0 
				% Check if structure start- and endpoints have individual offset uncertainties
				if length(SC.SIG.Support)>=ordPair(2) && ~isempty(SC.SIG.Support{ordPair(2)}) && isfield(SC.SIG.Support{ordPair(2)},[type{1} 'Offset'])
					error('Pitch angle errors can not be given explicitly if ''%s'' start and endpoints each have offset uncertainties.',type{1})
				end
				
				% Adjust support structure startpoint vertcial offset
				SC.RING{ordPair(1)}.([type{1} 'Offset'])(2) = SC.RING{ordPair(1)}.([type{1} 'Offset'])(2) - SC.RING{ordPair(1)}.([type{1} 'Roll'])(2)*structLength/2;
				% Adjust support structure endpoint vertcial offset
				SC.RING{ordPair(2)}.([type{1} 'Offset'])(2) = SC.RING{ordPair(2)}.([type{1} 'Offset'])(2) + SC.RING{ordPair(1)}.([type{1} 'Roll'])(2)*structLength/2;
			else
				% Get support structure pitch angle from vertcial start- and endpoint offsets
				SC.RING{ordPair(1)}.([type{1} 'Roll'])(2) = (SC.RING{ordPair(2)}.([type{1} 'Offset'])(2) - SC.RING{ordPair(1)}.([type{1} 'Offset'])(2))/structLength;
			end
			
			% Check if support structure yaw errors are given explicitly
			if SC.RING{ordPair(1)}.([type{1} 'Roll'])(3)~=0
				% Check if structure start- and endpoints have individual offset uncertainties
				if length(SC.SIG.Support)>=ordPair(2) && ~isempty(SC.SIG.Support{ordPair(2)}) && isfield(SC.SIG.Support{ordPair(2)},[type{1} 'Offset'])
					error('Yaw angle errors can not be given explicitly if ''%s'' start and endpoints each have offset uncertainties.',type{1})
				end
				
				% Adjust support structure startpoint horizontal offset
				SC.RING{ordPair(1)}.([type{1} 'Offset'])(1) = SC.RING{ordPair(1)}.([type{1} 'Offset'])(1) - SC.RING{ordPair(1)}.([type{1} 'Roll'])(3)*structLength/2;
				% Adjust support structure endpoint horizontal offset
				SC.RING{ordPair(2)}.([type{1} 'Offset'])(1) = SC.RING{ordPair(2)}.([type{1} 'Offset'])(1) + SC.RING{ordPair(1)}.([type{1} 'Roll'])(3)*structLength/2;
			else
				% Get support structure yaw angle from horizontal start- and endpoint offsets
				SC.RING{ordPair(1)}.([type{1} 'Roll'])(3) = (SC.RING{ordPair(2)}.([type{1} 'Offset'])(1) - SC.RING{ordPair(1)}.([type{1} 'Offset'])(1))/structLength;
			end
		end
	end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnet errors
function SC = applyMagnetError(SC,par)
	
	% Check if any uncertainties for the magnets are given
	if ~isfield(SC.SIG,'Mag')
		return;
	end
	% Loop over magnets
	for ord = SC.ORD.Magnet

		if isempty(SC.SIG.Mag{ord})
			continue
		end
		% Loop over uncertanties
		for field=fieldnames(SC.SIG.Mag{ord})'
			% Check if individual cutoff is given
			if iscell(SC.SIG.Mag{ord}.(field{1}))
				nSig = SC.SIG.Mag{ord}.(field{1}){2};
				sig = SC.SIG.Mag{ord}.(field{1}){1};
			else
				nSig = par.nSig;
				sig = SC.SIG.Mag{ord}.(field{1});
			end

			% Bending angle error gets applied differently
			if strcmp(field{1},'BendingAngle')
				SC.RING{ord}.BendingAngleError = sig * SCrandnc(nSig,1,1);
			else
				SC.RING{ord}.(field{1}) = sig .* SCrandnc(nSig,size(sig));
			end
		end
		
	end
end


