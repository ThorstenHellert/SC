function QuadOrds = getBPM2QuadPairing_ALSU_AR(SC)
% Gets the quadrupole-BPM pairing for BBA (for AR pretty simple)
			
	QuadOrds = SCgetOrds(SC.RING,'QF:|QD:|QFA:');
	
	QuadOrds = repmat(QuadOrds,2,1);

end
