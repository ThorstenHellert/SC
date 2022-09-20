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
%	Change of RF frequency in Hz
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'BPMords'` (`SC.ORD.BPM`)::
%    List of BPM ordinates at which the dispersion should be returned
% `'CAVords'` (`SC.ORD.CAV`)::
%    List of cavity ordinates with which the dispersion should be measured
% `'nSteps'` (2)::
%    Number of RF steps (1st RF step is considered the reference). If more than 2 steps are 
%    specified, the measurement is bi-directional
%
% RETURN VALUE
% ------------
% `eta`::
% 	Dispersion [m/Hz]
%
% EXAMPLES
% --------
% Calculate the dispersion with a 1 kHz rf frequency change.
% -------------------------------
% eta = SCgetDispersion(SC,1E3);
% -------------------------------
%
% SEE ALSO
% --------
% *SCsetCavs2SetPoints*



warning('TODO: needs to be updated!')


% Parse optional arguments
p = inputParser;
addOptional(p,'BPMords',SC.ORD.BPM);
addOptional(p,'CAVords',SC.ORD.Cavity);
addOptional(p,'nSteps',2);
parse(p,varargin{:});
par = p.Results;


% Define RF steps
for nCav=1:length(par.CAVords)
	RFsteps(nCav,:) = SC.RING{par.CAVords(nCav)}.FrequencySetPoint + linspace(-RFstep,RFstep,par.nSteps);
end


% Get reference BPM reading
Bref = reshape(SCgetBPMreading(SC,'BPMords',par.BPMords)',[],1);

if par.nSteps==2
	% Change RF frequency
	SC = SCsetCavs2SetPoints(SC,par.CAVords,'Frequency',RFstep,'add');
	
	% Calculate second BPM reading
	B = reshape(SCgetBPMreading(SC,'BPMords',par.BPMords)',[],1);
	
	% Calculate dispersion
	eta = (B-Bref)/RFstep;
else
	% Loop over frequency setpoints
	for nStep=1:par.nSteps
		
		% Change RF frequency
		SC = SCsetCavs2SetPoints(SC,par.CAVords,'Frequency',RFsteps(:,nStep),'abs');
		
		% Calculate BPM reading differences
		dB(nStep,:) = reshape(SCgetBPMreading(SC,'BPMords',par.BPMords)',[],1) - Bref;
	end
	
	% Linear regression
	eta = linspace(-RFstep,RFstep,par.nSteps)'\dB;
	
end



