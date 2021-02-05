function QuadOrds = getBPM2QuadPairing_ALSU_SR(SC,BPMords)
	
		
	% Define pairing
	BPM( 1).QuadName = 'QF1'; 
	BPM( 1).updown   = 'up'; 
	
	BPM( 2).QuadName = 'QD1'; 
	BPM( 2).updown   = 'up'; 
	
	BPM( 3).QuadName = 'QD1'; 
	BPM( 3).updown   = 'down'; 
	
	BPM( 4).QuadName = 'QF2'; 
	BPM( 4).updown   = 'up'; 	
	
	BPM( 5).QuadName = 'QF3'; 
	BPM( 5).updown   = 'up'; 
	
	BPM( 6).QuadName = 'QF4'; 
	BPM( 6).updown   = 'up'; 
	
	BPM( 7).QuadName = 'QF4'; 
	BPM( 7).updown   = 'down'; 
	
	BPM( 8).QuadName = 'QF5'; 
	BPM( 8).updown   = 'up'; 
	
	BPM( 9).QuadName = 'QF6'; 
	BPM( 9).updown   = 'up'; 

	BPM(10).QuadName = 'QF6'; 
	BPM(10).updown   = 'down'; 

	BPM(11).QuadName = 'QF6'; 
	BPM(11).updown   = 'up'; 

	BPM(12).QuadName = 'QF6'; 
	BPM(12).updown   = 'down'; 

	BPM(13).QuadName = 'QF5'; 
	BPM(13).updown   = 'up'; 

	BPM(14).QuadName = 'QF4'; 
	BPM(14).updown   = 'down'; 

	BPM(15).QuadName = 'QF2'; 
	BPM(15).updown   = 'up'; 

	BPM(16).QuadName = 'QF2'; 
	BPM(16).updown   = 'down'; 

	BPM(17).QuadName = 'QD1'; 
	BPM(17).updown   = 'up'; 

	BPM(18).QuadName = 'QF1'; 
	BPM(18).updown   = 'up'; 

	BPM(19).QuadName = 'QF1'; 
	BPM(19).updown   = 'down'; 
	
	for nDim=1:size(BPMords,1)
		i = 1;
		for ord=BPMords(nDim,:)
			BPMind  = str2num(SC.RING{ord}.FamName(4:end));
			tmpQuad = SCgetOrds(SC.RING,BPM(BPMind).QuadName);
			
			if isempty(tmpQuad)
				error('No quadrupole name could be found.')
			end
			
			switch BPM(BPMind).updown
				case 'up' % quad is upstream of BPM
					QuadOrds(nDim,i) = tmpQuad(find(tmpQuad>ord,1));
				case 'down' % quad is downstream of BPM
					QuadOrds(nDim,i) = tmpQuad(find(tmpQuad<ord,1,'last'));
			end
			i = i+1;
		end
	end
end
