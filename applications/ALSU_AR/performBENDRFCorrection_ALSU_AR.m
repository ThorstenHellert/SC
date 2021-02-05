function [SC,out] = performBENDRFCorrection_ALSU_AR(SC,targetF,DAsteps)
% Performs the linear optics correction required for a change in RF frequency for the ALS-U accumulator ring
	
	% Get estimated dipole setpoint change
	deltaF0  = SC.RING{SC.ORD.Cavity}.Frequency - targetF;
	deltaSP0 = 1-deltaF0*0.0019/1E3;
	
	% Set cavity and dipoles to new values
	SC = SCsetCavs2SetPoints(SC,SC.ORD.Cavity,'Frequency',-deltaF0,'add');
	SC = SCsetMags2SetPoints(SC,SCgetOrds(SC.RING,'BEND'),2,2,deltaSP0,'method','rel');
	SC = SCsetMags2SetPoints(SC,SCgetOrds(SC.RING,'QFA'),2,2,deltaSP0,'method','rel');
	
	% Get pre-LOCO lattice properties and save SC
	out.preBEND    = SCcalcLatticeProperties(SC,'DAsteps',DAsteps);
	out.preBEND.SC = SC;

	% Get modeal response matrix
	MCO = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'trackMode','ORB');
	
	% Define regularization parameters to try orbit feedback
	alphaVec = [10 5:-1:1];
	
	fprintf('Initial closed orbit deviation (hor//ver):   %.0fum // %.0fum\n',1E6*sqrt(mean(SCgetBPMreading(SC).^2,2)))
	
	% Decrease regularization parameter in orbit feedback loop
	for	alpha=alphaVec
		fprintf('Orbit correction with alpha = %d\n',alpha)
		
		% Get pseudo inverse
		MinvCO = SCgetPinv(MCO,'alpha',alpha);
		
		% Run feedback
		[CUR,ERROR] = SCfeedbackRun(SC,MinvCO,'target',0,'maxsteps',50,'verbose',1);
		if ERROR;warning('Feedback crashed!');break;end
		
		% Calculate intial and final rms BPM reading
		B0rms = sqrt(mean(SCgetBPMreading(SC ).^2,2));
		Brms  = sqrt(mean(SCgetBPMreading(CUR).^2,2));
		
		% Check if orbit feedback worked
		if mean(B0rms)<mean(Brms);break;end
		
		% Accept new machine
		SC = CUR;
		fprintf('CO mprovement with alpha = %d:\n hor: %.0fum -> %.0fum\n ver: %.0fum -> %.0fum\n',alpha,1E6*B0rms(1),1E6*Brms(1),1E6*B0rms(2),1E6*Brms(2))
	end
	
	% Perform LOCO
	SC = performLOCO_ALSU_AR(SC);
	
	% Get post-LOCO lattice properties and save SC
	out.postBEND    = SCcalcLatticeProperties(SC,'DAsteps',DAsteps);
	out.postBEND.SC = SC;
end
