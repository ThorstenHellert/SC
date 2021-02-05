function [SC,regNum] = performOrbitCorr_ALSU_SR(SC,RMstruct,varargin)
% performOrbitCorr_ALSU_SR
% ========================
%
% NAME
% ----
% performOrbitCorr_ALSU_SR - Performs orbit correction for the ALSU-SR
%
% SYNOPSIS
% --------
% `[SC,regNum] = performOrbitCorr_ALSU_SR(SC,RMstruct [,options])`
%
% DESCRIPTION
% -----------
% This function performs orbit correction in a loop with decreasing regularization parameter until
% no further reduction in BPM readings is found.
%
% INPUT
% -----
% `SC`::
% 	The SC base structure
% `RMstruct`::
% 	Structure containing all information relevant for orbit correction.
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'regVec'` (`RMstruct.alphaVec`)::
%	Regularization vector.
% `'regMode'` ('alpha')`)::
%	Regularization mode.
% `'checkTransmission'` (0)::
%	Check if for sufficient beam transmission is checked before each step.
% `'scaleDisp'` (`RMstruct.scaleDisp`)::
%	Scaling factor for dispersion 
% `'R0'` (zeros(2*length(RMstruct.BPMords),1))::
%	Reference orbit.
% `'BPMordsWeight'` (SCgetOrds(SC.RING,'BPM4|BPM5|BPM15|BPM16'))::
%	BPMs at which weighting factor should be applied.
% `'weight'` (10)::
%	Weighting factor for special BPMs.
% `'eps'` (1E-6)::
%	Noise level for feedback.
% `'verbose'` (1)::
%	If true, debug information is printed.
%
% RETURN VALUE
% ------------
% `SC`::
% 	The SC base structure with corrected orbit.
%
% SEE ALSO
% --------
% *performDiagCorr_ALSU_SR*, *performOrbitCorr_ALSU_SR*, *SCfeedbackRun*
	
	% Parse optional arguments
	p = inputParser;
	addOptional(p,'regVec',RMstruct.alphaVec);
	addOptional(p,'regMode','alpha');
	addOptional(p,'checkTransmission',0);
	addOptional(p,'scaleDisp',RMstruct.scaleDisp);
	addOptional(p,'R0',zeros(2*length(RMstruct.BPMords),1));
	addOptional(p,'weight',10);
	addOptional(p,'BPMordsWeight',SCgetOrds(SC.RING,'BPM4|BPM5|BPM15|BPM16'));
	addOptional(p,'eps',1E-6);
	addOptional(p,'verbose',1);
	parse(p,varargin{:});
	par = p.Results;
	
	% Define response matrix for orbit correction
	MCO = [RMstruct.RM RMstruct.scaleDisp*RMstruct.eta];
	

	if par.verbose;fprintf('Initial closed orbit deviation for all BPMs (hor//ver):   %.1fum // %.1fum\n',1E6*sqrt(mean(SCgetBPMreading(SC).^2,2)));end
	
	
	% Check if BPM weighting should be applied
	if isempty(par.weight)
		weight = ones(2,length(RMstruct.BPMords));
		weight = [weight(1,:)'; weight(2,:)'];
	else
		% Initialize weighting vecotr
		weight = ones(2,length(RMstruct.BPMords));
		
		% Get relevant indices of BPMs adjacent to chromatic sextupoles
		BPMordsSext = par.BPMordsWeight;
		ind = ismember(RMstruct.BPMords,BPMordsSext);
		% Set weight
		weight(:,ind) = par.weight;
		weight = [weight(1,:)'; weight(2,:)'];
	end
	
	% Apply weighting vector to RM
	MCO = MCO .* repmat(weight,1,size(MCO,2));
	
	% Decrease regularization parameter
	for	regNum=par.regVec
		
		if par.verbose;fprintf('Orbit correction with %s = %d\n',par.regMode,regNum);end
		
		% Get pseudo inverse
		MinvCO = SCgetPinv(MCO,par.regMode,regNum);
				
		% Check if closed orbit can be found
		SC = checkOrbit(SC);
		
		% Run feedback
		[CUR,ERROR] = SCfeedbackRun(SC,MinvCO,'weight',weight,'R0',par.R0,'target',0,'maxsteps',50,'scaleDisp',par.scaleDisp,'verbose',par.verbose,'BPMords',RMstruct.BPMords,'CMords',RMstruct.CMords,'eps',par.eps); 
		
		if ~ERROR
				
			if par.checkTransmission
				% Check for sufficient beam transmission for orbit mode
				maxTurns = SCgetBeamTransmission(CUR,...
					'nParticles',100,...
					'nTurns',100,...
					'verbose',par.verbose);
				
				% Check if transmission is below pseudo-orbit validity
				if maxTurns < 50
					fprintf('Transmission below 50 turns, abort orbit correction.\n')
					regNum = alphaVec(max([find(alphaVec==regNum)-1,1])); break
				end
			end
			
			% Get reference 
			Bref = reshape(par.R0,2,[]);
			% Get weighting
			tmpWeight = reshape(weight,[],2)';
			tmpWeight = tmpWeight / (norm(tmpWeight)/norm(ones(size(tmpWeight))));
			% Calculate initial and final BPM reading		
			B0rms = sqrt(mean(( (SCgetBPMreading(SC ,'BPMords',RMstruct.BPMords)-Bref) .* tmpWeight).^2,2));
			Brms  = sqrt(mean(( (SCgetBPMreading(CUR,'BPMords',RMstruct.BPMords)-Bref) .* tmpWeight).^2,2));
			
			% Check if orbit feedback worked
			if mean(B0rms) < mean(Brms)
				if par.verbose;fprintf('No further improvement with %s = %d\n',par.regMode,regNum);end
				regNum = par.regVec(max([find(par.regVec==regNum)-1,1])); break
			else
				SC = CUR;
				if par.verbose;fprintf('CO mprovement with %s = %d:\n hor: %.1fum -> %.1fum\n ver: %.1fum -> %.1fum\n',par.regMode,regNum,1E6*B0rms(1),1E6*Brms(1),1E6*B0rms(2),1E6*Brms(2));end
			end
		else
			regNum = par.regVec(max([find(par.regVec==regNum)-1,1])); break
		end
	end
	if par.verbose;fprintf('Final closed orbit deviation for all BPMs (hor//ver): %.1fum // %.1fum   with alpha = %d.\n',1E6*sqrt(mean(SCgetBPMreading(SC).^2,2)),regNum);end
	
end



% Check if closed orbit can be found
function SC = checkOrbit(SC)
	if  ~any(isnan(findorbit6(SC.RING,1))) 
		if ~strcmp(SC.INJ.trackMode,'ORB')
			SC.INJ.trackMode = 'ORB';
			fprintf('Switch to orbit mode\n')
		end
	else
		SC.INJ.trackMode  = 'pORB';
		SC.INJ.nTurns     = 100;
		SC.INJ.nParticles = 100;
		fprintf('Switch to pseudo orbit mode\n')
	end
end
