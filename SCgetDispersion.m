function eta = SCgetDispersion(SC,RFstep,varargin)
% SCgetDispersion
% ===============
%
% NAME
% ----
% SCgetDispersion - Measure dispersion by changing the rf frequency of the cavities
%
% SYNOPSIS
% --------
% `eta = SCgetDispersion(SC, RFstep)`
%
%
% DESCRIPTION
% -----------
% Calculates reference BPM reading, then changes the rf frequency and gets a
% second BPM reading in order to calculate the dispersion.
%
% INPUTS
% ------
% `SC`::
%	SC base structure
% `RFstep`::
%	Change of rf frequency
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'BPMords'` (`SC.ORD.BPM`):: List of BPM ordinates at which the dispersion should be returned
%
% RETURN VALUE
% ------------
% `eta`::
% 	Dispersion [m/Hz]
%
% EXAMPLES
% --------
% Calculate the dispersion with a 1kHz rf frequency change.
% -------------------------------
% eta = SCgetDispersion(SC,1E3);
% -------------------------------
%
% SEE ALSO
% --------
% *SCsetCavs2SetPoints*

% Parse optional arguments
p = inputParser;
addOptional(p,'BPMords',SC.ORD.BPM);
parse(p,varargin{:});

% Get reference BPM reading
Bref = reshape(SCgetBPMreading(SC,'BPMords',p.Results.BPMords)',[],1);

% Change RF frequency
SC = SCsetCavs2SetPoints(SC,SC.ORD.Cavity,'Frequency',RFstep,'add');

% Get second BPM reading
B = reshape(SCgetBPMreading(SC,'BPMords',p.Results.BPMords)',[],1);

% Calculate dispersion
eta = (B-Bref)/RFstep;
