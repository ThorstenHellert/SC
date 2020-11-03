function RING = SCscaleCircumference(RING,circ,varargin)
% SCscaleCircumference
% ====================
%
% Name
% ----
% SCscaleCircumference - scales the circumference of a ring
%
% Synopsis
% --------
% `RING = SCscaleCircumference(RING, circ [, mode])`
%
%
% Description
% -----------
% The circumference of `RING` is changed by scaling the lengths of the drift
% spaces in the lattice.
%
%
% MODE
% ----
% `'abs'` (default)::
%	The current circumference is scaled to the value `circ`.
% `'rel'`::
%	The current circumference is scaled by the value `circ`.
%
%
% Return Value
% ------------
% `RING`::
% 	Lattice with new circumference.

	% Get mode
	if nargin==3
		mode=varargin{1};
	else
		mode='abs';
	end

	% Length of the ring
	C=findspos(RING,length(RING)+1);

	% Sum of drift spaces in the ring
	D = 0;
	for ord=1:length(RING)
		if strcmp(RING{ord}.PassMethod,'DriftPass')
			D = D + RING{ord}.Length;
		end
	end
	switch mode
		case 'rel'
			% Scaling factor for the drifts
			Dscale = 1 - (1-circ) * C/D;
		case 'abs'
			Dscale = 1 - (C-circ)/D;
		otherwise
			error('Unsupported circumference scaling mode: ''%s''',mode)
	end

	% Apply scaling factor
	for ord=1:length(RING)
		if strcmp(RING{ord}.PassMethod,'DriftPass')
			RING{ord}.Length = RING{ord}.Length * Dscale;
		end
	end
end
