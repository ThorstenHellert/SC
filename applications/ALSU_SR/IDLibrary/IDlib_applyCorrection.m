function SC = IDlib_applyCorrection(SC,IDs,varargin)
% IDlib_applyCorrection
% =====================
%
% NAME
% ----
% IDlib_applyCorrection - Applys ID compensation
%
% SYNOPSIS
% --------
% `SC = IDlib_applyCorrection(SC,IDs [,options])`
%
% DESCRIPTION
% -----------
% This function closes the IDs and applies the quadrupole setpont change in order to
% compensate the linear optics pertubation by IDs and performes orbit-, tune- and chromaticity
% correction.
%
% INPUT
% -----
% `SC`::
% 	The SC base structure including IDs
% `IDs`::
% 	Cell string with ID names.
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'quadOrds'` (`SCgetOrds(SC.RING,'^QD1|^QF[1-2]|^QF[4-6]|BEND2|BEND3')`)::
%	Quadrupole ordinates for ID compensation.
% `'deltaK'` (`[]`)::
%	Quadrupole setpoint change.
% `'IDstrength'` (`ones(size(IDs))`)::
%	Scaling factor for IDs (only for series of SBENDs)
% `'IDmap'` (`''`)::
%	Name of kick map to be loaded.
% `'applyFeedbacks'` (`1`)::
%	Specify if orbit-, tune- and chromaticity correction is applied.
% `'matchChrom'` (`[2 1]`)::
%	Horizontal and vertical chromaticity target.
% `'chromOrds'` (`SCgetOrds(SC.RING,{'SF','SD'})`)::
%	Sextupole ordinates for chromaticity correction.
% `'matchTune'` (`[2 1]`)::
%	Horizontal and vertical tune target.
% `'tuneOrds'` (`SCgetOrds(SC.RING,{'^QF1','^QD1'})`)::
%	Quadrupole ordinates for tune correction.
% `'RMstruct'` (`[]`)::
%	Structure containing all relevant orbit correction parameters.
% `'verbose'` (`0`)::
%	If true debug information is printed.
%
% RETURN VALUE
% ------------
% `SC`::
% 	The SC base structure including closed and compensated IDs.
%
% SEE ALSO
% --------
% *IDLib_closeID*, *IDlib_includeIDs*, *IDlib_calcCorrection*, *performOrbitCorr_ALSU_SR*


	% Parse optional arguments
	p = inputParser;
	addOptional(p,'quadOrds',SCgetOrds(SC.RING,'^QD1|^QF[1-2]|^QF[4-6]|BEND2|BEND3'));
	addOptional(p,'deltaK',[]);
	addOptional(p,'IDstrength',ones(size(IDs)));
	addOptional(p,'IDmap',[]);
	addOptional(p,'applyFeedbacks',1);
	addOptional(p,'matchChrom',[2 1]);
	addOptional(p,'chromOrds',SCgetOrds(SC.RING,{'SF','SD'}));
	addOptional(p,'matchTune',[]);
	addOptional(p,'tuneOrds',SCgetOrds(SC.RING,{'^QF1','^QD1'}));
	addOptional(p,'RMstruct',[]); % RM used for orbit correction (inkl. dispersion), see: *performOrbitCorr_ALSU_SR*
	addOptional(p,'verbose',0);
	parse(p,varargin{:});
	par = p.Results;
	
	inputCheck()
	
	% Close IDs
	for nID=1:length(IDs)
		SC = IDLib_closeID(SC,IDs{nID},'IDstrength',par.IDstrength(nID),'IDmap',par.IDmap{nID},'PassMethod','ThinEPU2Pass');
	end
	
	% Apply quadrupole correction
	SC = SCsetMags2SetPoints(SC,par.quadOrds,2,2,par.deltaK,'method','add','dipCompensation',1);
	
	if par.applyFeedbacks
		% Calculate model response matrix
		if isempty(par.RMstruct)
			error('Response matrix structure for orbit feedback must be provided')
		end
		
		for n=1:3
			% Orbit feedback
			[SC,~] = performOrbitCorr_ALSU_SR(SC,par.RMstruct);
			
			% Match chromaticity
			SC = SClocoLib('fitChromaticity',SC,par.chromOrds,...
						   'targetChrom',par.matchChrom,...
						   'verbose',par.verbose);
			
			% Orbit feedback
			[SC,~] = performOrbitCorr_ALSU_SR(SC,par.RMstruct);
			
			% Tune feedback
			SC = SClocoLib('fitTune',SC,par.tuneOrds,...
				           'targetTune',par.matchTune,...
						   'TolX',1E-4,...
						   'TolFun',1E-3,...
						   'InitStepSize',[.02 .02],...
						   'verbose',par.verbose);
		end
		
		% Orbit feedback
		[SC,~] = performOrbitCorr_ALSU_SR(SC,par.RMstruct);
	end
	

	% Input check
	function inputCheck()
			
		if isempty(par.deltaK)
			error('No correction step ''deltaK'' specified.')
		end
		
		if length(par.quadOrds) ~= length(par.deltaK)
			error('Quadrupole ordinates must match size of ''deltaK''.')
		end
	end
end
