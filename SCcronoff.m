function RING = SCcronoff(RING,varargin)
% SCcronoff
% =========
%
% NAME
% ----
% SCcronoff - switch cavity/radiation on/off
%
% SYNOPSIS
% --------
% `RING = SCcronoff(RING, mode, ...)`
%
%
% DESCRIPTION
% -----------
% Depending on `mode` SCcronoff switches the cavities / the radiation in `RING`
% on or off. Possible `mode`s are `'radiationoff'`, `'radiationon'`, `'cavityoff'`,
% `'cavityon'`. Multiple modes can be specified.
%
%
% RETURN VALUE
% ------------
% `RING`::
% 	The modified base structure lattice.
%
% EXAMPLES
% --------
%
% Switch cavities and radiation in `SC.RING` off.
% ------------------------------------------------------------------
% SC.RING = SCcronoff(SC.RING,'cavityoff','radiatonoff');
% ------------------------------------------------------------------


	for i=1:(nargin-1)
		switch varargin{i}

			case 'radiationoff'
				for ord=1:length(RING)
					switch RING{ord}.PassMethod
						case 'BndMPoleSymplectic4RadPass'
							RING{ord}.PassMethod = 'BndMPoleSymplectic4Pass';
						case 'StrMPoleSymplectic4RadPass'
							RING{ord}.PassMethod = 'StrMPoleSymplectic4Pass';
					end
				end


			case 'radiationon'
				for ord=1:length(RING)
					switch RING{ord}.PassMethod
						case 'BndMPoleSymplectic4Pass'
							RING{ord}.PassMethod = 'BndMPoleSymplectic4RadPass';
						case 'StrMPoleSymplectic4Pass'
							RING{ord}.PassMethod = 'StrMPoleSymplectic4RadPass';
					end
				end


			case 'cavityoff'
				for ord=1:length(RING)
					if isfield(RING{ord},'Frequency')
						RING{ord}.PassMethod = 'IdentityPass';
					end
				end


			case 'cavityon'
				for ord=1:length(RING)
					if isfield(RING{ord},'Frequency')
						RING{ord}.PassMethod = 'RFCavityPass';
					end
				end


			otherwise
				warning('SCcronoff: mode %s not recognized. RING unchanged.',mode);
		end
	end
end
