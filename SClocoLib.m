function varargout = SClocoLib(funName,varargin)
% SClocoLib
% =========
%
% NAME
% ----
% SClocoLib - Function library to use LOCO with SC
%
% SYNOPSIS
% --------
% `varargout = SClocoLib(funName, varargin])`
%
%
% DESCRIPTION
% -----------
% `SClocoLib` is a function library intended to connect the workflows and data
% structures of the AT/MML function `loco` with the `SC` workflow. The input
% string `funName` defines the function which should be used. Additional input
% arguments and the return values of `SClocoLib` differ, depending on the
% called function. The following will describe each provided function in
% detail:
%
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% NAME
% ----
% setupLOCOmodel - Sets up the LOCO model
%
% SYNOPSIS
% --------
% `[LOCOmodel, LOCOflags, Init] = SClocoLib('setupLOCOmodel', SC [, options])`
%
% DESCRIPTION
% -----------
% This function sets up the LOCO model used by the function `loco` based on the
% `SC` structure.  Additionally lattice properties of `SC.RING` are calculated
% by *SCcalcLatticeProperties* and stored together with `SC` in the structure
% `Init`. This is usefull for further LOCO steps.
%
% INPUTS
% ------
% `SC`::
%	SC base structure.
%
% OPTIONS
% -------
% Additional arguments can be given as name-value pairs and are written as
% fields in `LOCOflags`.
%
% RETURN VALUE
% ------------
% `LOCOmodel`::
%	LOCO model used by `locoresponsematrix`.
% `LOCOflags`::
%	LOCO flags used by `loco`.
% `Init`::
%	Structure containing the initial `SC` structure and (disturbed) lattice
%	properties.
%
% EXAMPLES
% --------
% Set up the LOCO model, include dispersion and set the horizontal and vertical
% dispersion weights to 100.
% ------------------------------------------------------------------
% [LOCOmodel,LOCOflags,Init] = SClocoLib('setupLOCOmodel',SC,...
% 	'Dispersion','Yes',...
% 	'HorizontalDispersionWeight',.1E2,...
% 	'VerticalDispersionWeight',.1E2);
% ------------------------------------------------------------------
%
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% NAME
% ----
% getBPMCMstructure - Sets up the BPM and CM data structure for LOCO
%
% SYNOPSIS
% --------
% `[BPMData, CMData] =  SClocoLib('getBPMCMstructure', SC, CMstep [, options])`
%
% DESCRIPTION
% -----------
% This function sets up the BPM and CM data structure used by the function
% `loco` based on the `SC` registration.
%
% INPUTS
% ------
% `SC`::
%	SC base structure.
% `CMstep`::
%	CM step [rad] for `locoresponsematrix`.
%
% OPTIONS
% -------
% Additional arguments can be given as cell arrays of strings as
% type-name-value triples and are written as fields in in the BPM or CM structure
% (see examples) or to specify the used BPMs and CMs.
%
% RETURN VALUE
% ------------
% `BPMData`::
%	BPM data structure.
% `CMData`::
%	CM data structure.
%
% EXAMPLES
% --------
% Set up the BPM and CM data structures with a CM step of `0.1mrad` and include
% fitting the BPM gains and CM kicks and use only every second CM and BPM in both planes.
% ------------------------------------------------------------------
% [BPMData,CMData] =  SClocoLib('getBPMCMstructure',SC,1E-4,...
% 	{'BPM','FitGains','Yes'},...
% 	{'CM','FitKicks','Yes'},...
%   {'CMords',SC.ORD.CM{1}(1:2:end),SC.ORD.CM{2}(1:2:end)},...
%   {'BPMords',SC.ORD.BPM(1:2:end),SC.ORD.BPM(1:2:end)};
% ------------------------------------------------------------------
%
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% NAME
% ----
% getMeasurement - Simulates the response matrix measurement
%
% SYNOPSIS
% --------
% `LOCOmeasData =  SClocoLib('getMeasurement', SC, CMstep, RFstep [, options])`
%
% DESCRIPTION
% -----------
% Sets up the 'measured' response matrix data structure used by `loco`. Optional arguments
% are passed to *SCgetRespMat*.
%
% INPUTS
% ------
% `SC`::
%	SC base structure.
% `CMstep`::
%	CM step [rad] for *SCgetRespMat*.
% `RFstep`::
%	RF frequency step [Hz] for *SCgetDispersion*.
%
% RETURN VALUE
% ------------
% `LOCOmeasData`::
%	Data structure containing the response matrix.
%
% EXAMPLES
% --------
% Get the orbit response matrix and the dispersion measurement using CM steps
% of 0.1mrad and an rf frequency step of 1kHz, respectively.
% ------------------------------------------------------------------
% LOCOmeasData =  SClocoLib('getMeasurement',SC,1E-4,1E3);
% ------------------------------------------------------------------
% Get the orbit response matrix and the dispersion measurement using a BPM variation
% of 0.1mm and an rf frequency step of 1kHz, respectively and save the used CM steps.
% ------------------------------------------------------------------
% [LOCOmeasData, CMsteps] =  SClocoLib('getMeasurement',SC,1E-4,1E3,'mode','fixedOffset');
% ------------------------------------------------------------------
%
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% NAME
% ----
% setupFitparameters - Set up the LOCO fit parameter structure
%
% SYNOPSIS
% --------
% `FitParameters = SClocoLib('setupFitparameters', SC, RING0, LOCOmodel, RFstep [, parameters])`
%
% DESCRIPTION
% -----------
% Sets up the fit parameter structure which shall be used by `loco`.
% Additionally the inital values of the fit parameters are stored, which is
% needed in order to eventually apply the lattice correction.
%
% INPUTS
% ------
% `SC`::
%	SC base structure.
% `RING0`::
%	Initial (disturbed) lattice cell structure which was used for the response matrix measurement.
% `LOCOmodel`::
%	LOCO model structure.
% `RFstep`::
%	RF frequency step [Hz] for `locoresponsematrix`.
%
% PARAMETERS
% ----------
% Additional arguments are given as cell arrays and specify the
% fit parameters.  Each cell must be given as {`ordinates`, `normal/skew`,
% `individual/family`, `deltaK`} quadrupel, see examples.
%
% RETURN VALUE
% ------------
% `FitParameters`::
%	Fit parameter structure for `loco`.
%
% EXAMPLES
% --------
% Set up the LOCO fit parameter structure using all `QF` and `QD` normal
% quadrupoles (see also *SCgetOrds*) which are individually powered and using a
% strength variation of 1E-3 and 1E-4 to calculate the derivatives,
% respectively. The disturbed lattice `RING0` is used to identify the inital
% setpoints of the fit parameters and an rf frequency step of 1kHz is assumed.
% ------------------------------------------------------------------
% ordQF = SCgetOrds(SC.RING,'QF');
% ordQD = SCgetOrds(SC.RING,'QD');
% FitParameters = SClocoLib('setupFitparameters',SC,RING0,LOCOmodel,1E3,...
% 	{ordQF,'normal','individual',1E-3},...
% 	{ordQD,'normal','individual',1E-4});
% ------------------------------------------------------------------
%
% Set up the LOCO fit parameter structure using all `QFA` normal quadrupoles
% which are powered as a group (one fit parameter), using a strength variation
% of 1E-4 to calculate the derivatives. The previously described structure
% `Init` (see `'setupLOCOmodel'`) is used to identify the inital setpoints of
% the fit parameters and an rf frequency step of 1kHz is assumed.
% ------------------------------------------------------------------
% FitParameters = SClocoLib('setupFitparameters',SC,Init.SC.RING,LOCOmodel,1E3,...
% 	{SCgetOrds(SC.RING,'QFA'),'normal','family',1E-4},...
% ------------------------------------------------------------------
%
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
%
% NAME
% ----
% applyLatticeCorrection - Applies the LOCO lattice correction
%
% SYNOPSIS
% --------
% `SC = SClocoLib('applyLatticeCorrection', SC, FitParameters)`
%
% DESCRIPTION
% -----------
% Applies the calculated lattice correction by adjusting the setpoints in
% `SC.RING`.
%
% INPUTS
% ------
% `SC`::
%	SC base structure.
% `FitParameters`::
%	Fit parameter structure (calculated by `loco`).
%
% RETURN VALUE
% ------------
% `SC`::
%	SC base structure with applied lattice correction.
%
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
%
% NAME
% ----
% applyDiagnosticCorrection - Applies the LOCO result to BPMs and CMs
%
% SYNOPSIS
% --------
% `SC = SClocoLib('applyDiagnosticCorrection', SC, CMstep, CMData, BPMData [, options])`
%
% DESCRIPTION
% -----------
% Corrects the BPM and CM errors in `SC` by applying the LOCO fit results. If quadrupoles are used
% as CMs, the actual CM kick differs from the setpoint due to feed down effects. The (linear)
% routine used by LOCO to calculate the response matrices does not include that effect. Thus, the
% fitted CM calibration will differ from the actual calibration error of the CMs. This effect can be
% mitigated by using the optional 'CMcalOffsets'. The option 'meanToZero' helps to deal with cases
% where e.g. all fitted horizontal CM kicks are 2% too large while all horizontal BPM gains are 2%
% too low.
% 
% INPUTS
% ------
% `SC`::
%	SC base structure.
% `CMstep`::
%	CM step [rad] used in `getBPMCMstructure`.
% `BPMData`::
%	BPM data structure.
% `CMData`::
%	CM data structure.
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'CMcalOffsets'` ([])::
%   1x2 cell aray containg the horizontal and vertical CM calibration offsets
% `'meanToZero'` (0)::
%   If true, the fitted calibration factors are subtracted by their mean value
%   before the correction is aplied
% `'outlierRemovalAt'` ([])::
%	If this option is set, any fitted calibration error above the threshold is
%   discarded and not applied to the SC structure.
%
% RETURN VALUE
% ------------
% `SC`::
%	SC base structure with applied correction to BPM and CM errors.
%
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
%
% NAME
% ----
% applyOrbitCorrection - Applies orbit correction
%
% SYNOPSIS
% --------
% `SC = SClocoLib('applyOrbitCorrection', SC [, options])`
% 
% DESCRIPTION
% -----------
% Applies orbit feedback using *SCfeedbackRun*.
%
% INPUTS
% ------
% `SC`::
%	SC base structure.
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'Minv'` ([])::
%	Pseudo inverse of the response matrix. If not given explicitly,
%	*SCgetModelRM* is called and the Tikhonov regularization is used by
%   *SCgetPinv* to calculate `Minv`.
% `'alpha'` (50)::
%   Tikhonov regularization parameter.
% `'CMords'` (SC.ORD.CM)::
%   CM oridnates used for orbit correction.
% `'BPMords'` (SC.ORD.BPM)::
%   BPM oridnates used for orbit correction.
% `'verbose'` (0)::
%	If true, debug information is printed.
%
% RETURN VALUE
% ------------
% `SC`::
%	SC base structure with applied orbit correction.
%
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% NAME
% ----
% fitChromaticity - Fit chromaticity
%
% SYNOPSIS
% --------
% `SC = SClocoLib('fitChromaticity', SC, sFamOrds, [, options])`
%
% DESCRIPTION
% -----------
% Applies a chromaticity correction using two sextupole families. The absolute initial 
% setpoint variation wihtin one family remains unchanged. Note:
% this is not beam based but assumes the chromaticities can be measured reasonably well.
%
% INPUTS
% ------
% `SC`::
%	SC base structure.
% `sFamOrds`::
%	[1x2] cell array of sextupole ordinates used for matching
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'targetChrom'` ([])::
%   [1x2] array of target chromaticities. If not specified, the chromaticities as
%   calcualted from `SC.IDEALRING` are being used.
% `'InitStepSize'` ([2 2])::
%   Initial step size for `fminsearch`
% `'TolX'` (1E-4)::
%   Step size tolerance used by `fminsearch`
% `'TolFun'` (1E-3)::
%   Merrit function tolerance used by `fminsearch`
% `verbose` (0)::
%	If true, debug information is printed.
%
% RETURN VALUE
% ------------
% `SC`::
%	SC base structure with corrected chromaticity.
%
% EXAMPLES
% --------
% Match the chromaticity in both planes to 1 using all magnets named `'SF'` and
% `'SD'` (see also *SCgetOrds*).
% ------------------------------------------------------------------
% SC = SClocoLib('fitChromaticity',SC,SCgetOrds(SC.RING,{'SF','SD'}),[1 1]);
% ------------------------------------------------------------------
%
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% NAME
% ----
% fitTune - Fit tunes
%
% SYNOPSIS
% --------
% `SC = SClocoLib('fitTune', SC, qFamOrds, [, options])`
%
% DESCRIPTION
% -----------
% Applies a tune correction using two quadrupole families. The absolute initial 
% setpoint variation wihtin one family, e.g. from LOCO remains unchanged. Note:
% this is not beam based but assumes the tunes can be measured reasonably well.
%
% INPUTS
% ------
% `SC`::
%	SC base structure.
% `qFamOrds`::
%	[1x2] cell array of quadrupole ordinates used for matching
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'FitInteger'` (1)::
%   Flag specifying if the integer part should be fitted as well.
% `'targetTune'` ([])::
%   [1x2] array of target tunes. If not specified, the tunes as
%   calcualted from `SC.IDEALRING` are being used.
% `'TolX'` (1E-4)::
%   Step size tolerance used by `fminsearch`
% `'InitStepSize'` ([.01 .01])::
%   Initial step size for `fminsearch`
% `'TolFun'` (1E-3)::
%   Merrit function tolerance used by `fminsearch`
% `verbose` (0)::
%	If true, debug information is printed.
%
% RETURN VALUE
% ------------
% `SC`::
%	SC base structure with corrected chromaticity.
%
% EXAMPLES
% --------
% Match the fractional tunes in both planes to the ideal tunes using all 
% magnets named `'QF'` and `'QD'` (see also *SCgetOrds*).
% ------------------------------------------------------------------
% SC = SClocoLib('fitTune',SC,SCgetOrds(SC.RING,{'QF','QD'}),'FitInteger',0);
% ------------------------------------------------------------------
%
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
% NAME
% ----
% plotStatus - Plots the current LOCO correction
%
% SYNOPSIS
% --------
% `SClocoLib('plotStatus', SC, Init, BPMData, CMData)`
%
% DESCRIPTION
% -----------
% Plots the current beta beat, dispersion error and BPM and CM calibration and
% roll errors.
%
% INPUTS
% ------
% `SC`::
%	SC base structure.
% `Init`::
%	Structure containing the initial `SC` structure and (disturbed) lattice
%	properties.
% `BPMData`::
%	BPM data structure.
% `CMData`::
%	CM data structure.


	% Call function
	eval([funName,'(varargin{:})']);

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Set up LOCO model
	function setupLOCOmodel(SC,varargin)

		% Create LOCO model from ideal ring
		LOCOmodel.CavityFrequency   = SC.IDEALRING{SC.ORD.Cavity}.Frequency;
		LOCOmodel.CavityHarmNumber  = SC.IDEALRING{SC.ORD.Cavity}.HarmNumber;
		LOCOmodel.Lattice           = SC.IDEALRING;
		
		% Calculate disturbed lattice properties and save initial machine state
		Init.SC      = SC;

		% Define LOCO parameters
		LOCOflags.HorizontalDispersionWeight    = 10;    % Hor. dispersion VS. ORM elements
		LOCOflags.VerticalDispersionWeight      = 10;    % Ver. dispersion VS. ORM elements
		LOCOflags.Dispersion                    = 'No';  % Include dispersion
		LOCOflags.FitHCMEnergyShift             = 'No';  % Fit HCM energy Shift
		LOCOflags.FitVCMEnergyShift             = 'No';  % Fit VCM energy Shift
		LOCOflags.SVmethod                      = 1E-3;  % Cut off
		LOCOflags.AutoCorrectDelta              = 'No';  % NO!
		LOCOflags.Normalize                     = 'Yes'; % Normalization flag
		LOCOflags.Linear                        = 'Yes'; % Response matrix calculator
		LOCOflags.Coupling	                    = 'No'; % Include off-diagonal ORM elements
		LOCOflags.Dispersion                    = 'No'; % Include dispersion

		% Set name/pair-values in LOCO flags structure
		for i=1:2:(length(varargin)-1)
			LOCOflags.(varargin{i}) = varargin{i+1};
		end

		% Define output arguments
		varargout{1} = LOCOmodel;
		varargout{2} = LOCOflags;
		varargout{3} = Init;

	end

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Set up BPM and CM data structure structures
	function getBPMCMstructure(SC,CMsteps,varargin)

		% Default BPM data
		BPMData.FamName           = 'BPM';
		BPMData.BPMIndex          = SC.ORD.BPM(:);
		BPMData.HBPMIndex         = [1:length(SC.ORD.BPM)]';
		BPMData.VBPMIndex         = [1:length(SC.ORD.BPM)]';
		BPMData.FitGains          = 'No';

		% Default CM data
		CMData.FamName          = 'CM';
		CMData.HCMIndex         = SC.ORD.CM{1}(:);
		CMData.VCMIndex         = SC.ORD.CM{2}(:);
		CMData.FitKicks         = 'No';
		
		% Set name/pair-values in BPM/CM data structure
		for i=1:(length(varargin))
			switch varargin{i}{1}
				case 'CMords'
					CMData.HCMIndex     = varargin{i}{2}(:);
					CMData.VCMIndex     = varargin{i}{3}(:);
				case 'BPMords'
					BPMData.BPMIndex  = unique([varargin{i}{2}(:),varargin{i}{3}(:)]);
					BPMData.HBPMIndex = find(ismember(BPMData.BPMIndex,varargin{i}{2}));
					BPMData.VBPMIndex = find(ismember(BPMData.BPMIndex,varargin{i}{3}));
				case 'BPM'
					BPMData.(varargin{i}{2}) = varargin{i}{3};
				case 'CM'
					CMData.(varargin{i}{2}) = varargin{i}{3};
				otherwise
					error('Unsuported type.')
			end
		end
		
		
		% If CM step is given as single number, expand to the right dimensions
		if isnumeric(CMsteps) && length(CMsteps)==1
			CMsteps = {repmat(CMsteps,length(CMData.HCMIndex),1),repmat(CMsteps,length(CMData.VCMIndex),1)};
		end
		CMData.CMsteps = CMsteps;

		% Set other default values
		BPMData.HBPMGain      = ones(size(BPMData.HBPMIndex));
		BPMData.VBPMGain      = ones(size(BPMData.VBPMIndex));
		BPMData.HBPMCoupling  = zeros(size(BPMData.HBPMIndex));
		BPMData.VBPMCoupling  = zeros(size(BPMData.VBPMIndex));
		CMData.HCMCoupling    = zeros(size(CMData.HCMIndex));
		CMData.VCMCoupling    = zeros(size(CMData.VCMIndex));
		
		
		% Set CM kicks in LOCO structure
		CMData.HCMKicks = 1E3*2*CMData.CMsteps{1}(:).*ones(size(CMData.HCMIndex)); % [mrad]
		CMData.VCMKicks = 1E3*2*CMData.CMsteps{2}(:).*ones(size(CMData.VCMIndex)); % [mrad]

		% Define output arguments
		varargout{1} = BPMData;
		varargout{2} = CMData;
	end

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get orbit response matrix and dispersion from measurement
	function getMeasurement(SC,CMstep,deltaRF,BPMords,CMords,varargin)

		LocoMeasData.RF         = SC.RING{SC.ORD.Cavity(1)}.Frequency;
		LocoMeasData.DeltaRF    = deltaRF;
		LocoMeasData.BPMSTD     = 1E-3*ones(1,2*length(BPMords)); % [mm]
				
		
		[RM,Err,CMsteps] = SCgetRespMat(SC,CMstep,BPMords,CMords,varargin{:});
		
		% Take only absolute amplitude and transpose to match LOCO data
		CMsteps = cellfun(@transpose,cellfun(@max,cellfun(@abs,CMsteps,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false);
		
		LocoMeasData.M      = 2*1000* repmat([CMsteps{1};CMsteps{2}]',size(RM,1),1) .* RM;
% 		LocoMeasData.BPMSTD = Err;
		
		LocoMeasData.Eta = LocoMeasData.DeltaRF*1000 * SCgetDispersion(SC,LocoMeasData.DeltaRF,'BPMords',BPMords,'nSteps',3);

		% Define output arguments
		varargout{1} = LocoMeasData;
		varargout{2} = CMsteps;

	end

	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Set up LOCO fit parameter structure
	function setupFitparameters(SC,RING0,MODEL,DeltaRF,varargin)

		FitParameters.FitRFFrequency = 'Yes';
		FitParameters.DeltaRF        = DeltaRF;

		nGroup=1;
		for nFam = 1:length(varargin)
			% Reset element index
			nElem = 1;

			% Loop over family ordinates
			for ord = varargin{nFam}{1}
				% Swith type
				switch varargin{nFam}{2}
					case 'normal'
						FitParameters.Params{nGroup,1}(nElem).FieldName    = 'PolynomB';
						FitParameters.Params{nGroup,1}(nElem).SCFieldName  = 'SetPointB';
					case 'skew'
						FitParameters.Params{nGroup,1}(nElem).FieldName   = 'PolynomA';
						FitParameters.Params{nGroup,1}(nElem).SCFieldName = 'SetPointA';
					otherwise
						error('Unsoported type.')
				end
				% Create structure for LOCO
				FitParameters.Params{nGroup,1}(nElem).ElemIndex  = ord;
				FitParameters.Params{nGroup,1}(nElem).FieldIndex = {1,2};
				FitParameters.Params{nGroup,1}(nElem).Function   = @(x) x;
				FitParameters.Params{nGroup,1}(nElem).Args       = {};
				% Get fit parameter values
				FitParameters.Values(nGroup,1)      = MODEL.Lattice{ord}.(FitParameters.Params{nGroup}(nElem).FieldName)(2);
				FitParameters.IdealValues(nGroup,1) = SC.IDEALRING{ ord}.(FitParameters.Params{nGroup}(nElem).FieldName)(2);
				FitParameters.OrigValues(nGroup,1)  = RING0{ord}.(        FitParameters.Params{nGroup}(nElem).SCFieldName)(2);

				% Define initial deltas
				FitParameters.Deltas(nGroup,1) = varargin{nFam}{4};

				% Determine if group or element index gets increased
				if strcmp(varargin{nFam}{3},'family')
					nElem = nElem + 1;
					if ord==varargin{nFam}{1}(end)
						nGroup = nGroup + 1;
					end
				else
					nGroup = nGroup + 1;
				end
			end
		end

		% Define output arguments
		varargout{1} = FitParameters;
	end

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Apply lattice correction step
	function applyLatticeCorrection(SC,FitParameters)
		% Loop over fit parameters
		for nGroup = 1:length(FitParameters.Params)
			% Loop over elements in group
			for nElem = 1:length(FitParameters.Params{nGroup})
				ord      = FitParameters.Params{nGroup}(nElem).ElemIndex;
				field    = FitParameters.Params{nGroup}(nElem).SCFieldName;
				setpoint = FitParameters.OrigValues(nGroup) + (FitParameters.IdealValues(nGroup)-FitParameters.Values(nGroup));
			
				% Apply setpoint correction
				switch field
					case 'SetPointB' % Normal quadrupole
						SC = SCsetMags2SetPoints(SC,ord,2,2,setpoint,'dipCompensation',1);
					case 'SetPointA' % Skew quadrupole
						SC = SCsetMags2SetPoints(SC,ord,1,2,setpoint);
				end
			end
		end
		% Update magnets
		SC = SCupdateMagnets(SC,SC.ORD.Magnet);

		% Define output arguments
		varargout{1} = SC;
	end

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Apply BPM and CM corrections
	function applyDiagnosticCorrection(SC,CMstep,CMData,BPMData,varargin)
		% Parse optional arguments
		p = inputParser;
		addOptional(p,'CMcalOffsets',[]);
		addOptional(p,'meanToZero',0);
		addOptional(p,'outlierRemovalAt',[]);
		parse(p,varargin{:});
				
		if strcmp(CMData.FitCoupling,'Yes')
			warning('CM roll correction not implemented yet.')
		end
		
		% Check if CM step is specified for each CM
		if (~iscell(CMstep) && ~length(CMstep)==1 ) || (iscell(CMstep) && ( length(CMstep{1})~=length(CMData.HCMIndex) && length(CMstep{2})~=length(CMData.VCMIndex) ))
			error('CM steps must be defined as single value or cell array matching the number of used HCM and VCM.')
		end
		
		% Expand CM steps to cell array
		if ~iscell(CMstep)
			CMstep = {repmat(CMstep,length(CMData.HCMIndex),1),repmat(CMstep,length(CMData.VCMIndex),1)};
		end
				
		% Apply CM calibration error correction
		fields={'H','V'};SCfields={'CalErrorB','CalErrorA'};
		for nDim=1:2
			if strcmp(CMData.FitKicks,'Yes')
				fitCalCM{ nDim} = CMData.([fields{nDim} 'CMKicks'])./CMstep{nDim}/1000/2;

				% Add CM calibration offset
				if ~isempty(p.Results.CMcalOffsets)
					fitCalCM{nDim} = fitCalCM{nDim} - p.Results.CMcalOffsets{nDim};
				end
				
				% Check for outlier removal
				if ~isempty(p.Results.outlierRemovalAt)
					fitCalCM{nDim}(abs(1-fitCalCM{nDim})>=p.Results.outlierRemovalAt) = 1;
				end
				
				% Set mean CM calibration to zero
				if p.Results.meanToZero==1
					fitCalCM{nDim} = fitCalCM{nDim} + mean(1-fitCalCM{nDim});
				end
				
				% Apply CM calibration correction
				i=1;
				for ord=CMData.([fields{nDim} 'CMIndex'])'
					SC.RING{ord}.(SCfields{nDim})(1) = SC.RING{ord}.(SCfields{nDim})(1) + (1-fitCalCM{nDim}(i));
					i=i+1;
				end
			end
			
			% Apply BPM calibration error correction
			if strcmp(BPMData.FitGains,'Yes')
				% Get fitted BPM calibration errors
				fitCalBPM  = BPMData.([fields{nDim} 'BPMGain']);

				% Check for outlier removal
				if ~isempty(p.Results.outlierRemovalAt)
					fitCalBPM(abs(1-fitCalBPM)>=p.Results.outlierRemovalAt) = 1;
				end

				% Set mean BPM calibration to zero
				if p.Results.meanToZero==1
					fitCalBPM = fitCalBPM + mean(1-fitCalBPM);
				end
							
				% Apply BPM calibration error correction
				i=1;
				for ord=BPMData.BPMIndex(BPMData.([fields{nDim} 'BPMIndex']))'
					SC.RING{ord}.CalError(nDim) = SC.RING{ord}.CalError(nDim) + (1-fitCalBPM(i));
					i=i+1;
				end
			end
		end
				
		% Apply BPM roll error correction
		if strcmp(BPMData.FitCoupling,'Yes')
			% Get fitted BPM rolls
			fitRollBPM = [BPMData.VBPMCoupling-BPMData.HBPMCoupling]'/2;
			
			% Apply BPM calibration error correction
			i=1;
			for ord=BPMData.BPMIndex'
				SC.RING{ord}.Roll = SC.RING{ord}.Roll - fitRollBPM(i);
				i=i+1;
			end
		end
		% Define output arguments
		varargout{1} = SC;
	end

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Perform orbit correction
	function applyOrbitCorrection(SC,varargin)
		% Parse optional arguments
		p = inputParser;
		addOptional(p,'Minv',[]);
		addOptional(p,'alpha',50);
		addOptional(p,'CMords',SC.ORD.CM);
		addOptional(p,'BPMords',SC.ORD.BPM);
		addOptional(p,'verbose',0);
		parse(p,varargin{:});
		par=p.Results;

		if isempty(par.Minv)
			% Calculate orbit reponse matrix
			M = SCgetModelRM(SC,par.BPMords,par.CMords,'trackMode','ORB');
			% Check if something went wrong
			if any(isnan(M(:)));error('NaN in model response, aborting.');end
			% Get pseudo inverse
			par.Minv = SCgetPinv(M,'alpha',par.alpha);
		end
		% Apply orbit feedback
		[CUR,ERROR] = SCfeedbackRun(SC,par.Minv,'target',0,'maxsteps',30,'BPMords',par.BPMords,'CMords',par.CMords,'verbose',par.verbose);
		if ERROR;warning('Feedback crashed.');else;SC=CUR;end

		% Define output arguments
		varargout{1} = SC;
	end

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Match chromaticity (so far not beam based)
	function matchChromaticity(SC,sFamOrds,chromTarget)
		
		warning('This function will be removed in a future release. Use ''fitChromaticity'' within this library instead')
		
		% Match chromaticity
		[tmp,~]=atmatchchromdelta(SC.RING,chromTarget(:),sFamOrds);

		% Apply fit to setpoints
		for nFam=1:length(sFamOrds)
			factor = tmp{sFamOrds{nFam}(1)}.PolynomB(3)/SC.RING{sFamOrds{nFam}(1)}.PolynomB(3);
			for ord=sFamOrds{nFam}
				SC.RING{ ord}.SetPointB(3) = SC.RING{ ord}.SetPointB(3) * factor;
			end
		end
		% Update magnets
		SC  = SCupdateMagnets(SC,SC.ORD.Magnet);

		% Define output arguments
		varargout{1} = SC;

	end

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Fit chromaticity (not beam based)
	function fitChromaticity(SC,sOrds,varargin)
		
		% Parse optional arguments
		p = inputParser;
		addOptional(p,'targetChrom',[]);
		addOptional(p,'verbose',0);
		addOptional(p,'InitStepSize',[2 2]);
		addOptional(p,'TolX',1E-4);
		addOptional(p,'TolFun',1E-3);
		addOptional(p,'sepTunesWithOrds',[]);
		addOptional(p,'sepTunesDeltaK',[]);
		parse(p,varargin{:});
		par = p.Results;
		
		% Select target chromaticity
		if isempty(par.targetChrom)
			[~, ~, par.targetChrom] = atlinopt(SC.IDEALRING,0,[]);
		end
		
		% Check if something wnet wrong
		if any(isnan(par.targetChrom))
			fprintf('Target chromaticity must not contain NaN. Aborting.\n')
			varargout{1} = SC;
			return
		end
		
		% Copy intial SC state (for convinience)
		SC0 = SC;
		
		% Seperate tunes for fitting
		if ~isempty(par.sepTunesWithOrds) && ~isempty(par.sepTunesDeltaK)
			for nFam=1:length(par.sepTunesWithOrds)
				SC = SCsetMags2SetPoints(SC,par.sepTunesWithOrds{nFam},2,2,par.sepTunesDeltaK(nFam),'method','add');
			end
		end
		
		% Print initial conditions
		if par.verbose
			[~, ~, xi0] = atlinopt(SC.RING,0,[]); fprintf('Fitting chromaticities from [%.3f,%.3f] to [%.3f,%.3f].\n',xi0,par.targetChrom)
		end
		
		% Initial setpoints
		for nFam=1:length(sOrds)
			for n=1:length(sOrds{nFam})
				SP0{nFam}(n) = SC.RING{sOrds{nFam}(n)}.SetPointB(3);
			end
		end
		
		% Define fit function
		fun = @(x) fitFunction(SC,sOrds,x,SP0,par.targetChrom);
		
		% Run solver
		sol = fminsearch(fun,par.InitStepSize,optimset('TolX',par.TolX,'TolFun',par.TolFun));
		
		% Apply solution
		SC = applySetpoints(SC0,sOrds,sol,SP0);
		
		% Print final falues
		if par.verbose
			[~, ~, xi1] = atlinopt(SC.RING,0,[]); fprintf('        Final chromaticity: [%.3f,%.3f]\n          Setpoints change: [%.2f,%.2f]\n',xi1,sol)
		end
		
		% Define output arguments
		varargout{1} = SC;

		% Define the fitting function
		function out = fitFunction(SC,qOrds,setpoints,SP0,target)
			% Apply setpoints
			SC = applySetpoints(SC,qOrds,setpoints,SP0);
			
			% Get chromaticity
			[~, ~, xi] = atlinopt(SC.RING,0,[]);
			
			% Get figure of merrit
			out = sqrt(mean( (xi(:) - target(:)).^2 ));
		end
		
		% Apply setpoints
		function SC = applySetpoints(SC,ords,setpoints,SP0)	
			for nF=1:length(ords)
				SC = SCsetMags2SetPoints(SC,ords{nF},2,3,setpoints(nF)+SP0{nF},'method','abs');
			end
		end
		
	end

	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot status of the correction chain
	function fitTune(SC,qOrds,varargin)
		
		% Parse optional arguments
		p = inputParser;
		addOptional(p,'targetTune',[]);
		addOptional(p,'verbose',0);
		addOptional(p,'TolX',1E-4);
		addOptional(p,'TolFun',1E-3);
		addOptional(p,'InitStepSize',[.01 .01]);
		addOptional(p,'FitInteger',1);
		parse(p,varargin{:});
		par = p.Results;
		
		% Select target tune
		if isempty(par.targetTune)
			if par.FitInteger
				[ld, ~, ~]  = atlinopt(SC.IDEALRING,0,1:length(SC.IDEALRING)+1);
				par.targetTune = ld(end).mu/2/pi;
			else
				[~, par.targetTune, ~]  = atlinopt(SC.IDEALRING,0);
			end
		end
		
		% Print initial conditions
		if p.Results.verbose
			nu0 = getLatProps(SC,par.FitInteger); fprintf('Fitting tunes from [%.4f,%.4f] to [%.4f,%.4f].\n',nu0,par.targetTune)
		end
		
		% Initial setpoints
		for nFam=1:length(qOrds)
			for n=1:length(qOrds{nFam})
				SP0{nFam}(n) = SC.RING{qOrds{nFam}(n)}.SetPointB(2);
			end
		end

		% Define fit function
		fun = @(x) fitFunction(SC,qOrds,x,SP0,par.targetTune,par.FitInteger);
		
		% Run solver
		sol = fminsearch(fun,par.InitStepSize,optimset('TolX',par.TolX,'TolFun',par.TolFun));
		
		% Apply solution
		SC = applySetpoints(SC,qOrds,sol,SP0);
		
		% Print final falues
		if p.Results.verbose
			nu1 = getLatProps(SC,par.FitInteger); fprintf('       Final tune: [%.4f,%.4f]\n  Setpoints change: [%.3f,%.3f]\n',nu1,sol)
		end
		
		% Define output arguments
		varargout{1} = SC;

		% Define the fitting function
		function out = fitFunction(SC,qOrds,setpoints,SP0,target,FitInteger)
			% Apply setpoints
			SC = applySetpoints(SC,qOrds,setpoints,SP0);
			
			% Calcualte lattice properties
			nu = getLatProps(SC,FitInteger);
			
			% Get figure of merrit
			out = sqrt(mean( (nu(:) - target(:)).^2 ));
		end
		
		% Get tunes
		function nu = getLatProps(SC,FitInteger)
			if FitInteger
				[tmp, ~, ~]  = atlinopt(SC.RING,0,1:length(SC.RING)+1);
				nu = tmp(end).mu/2/pi;
			else
				[~, nu, ~]  = atlinopt(SC.RING,0);
			end
		end
		
		% Apply setpoints
		function SC = applySetpoints(SC,ords,setpoints,SP0)
			for nFam=1:length(ords)
				SC = SCsetMags2SetPoints(SC,ords{nFam},2,2,setpoints(nFam)+SP0{nFam},'method','abs','dipCompensation',1);
			end
		end

	end



	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot status of the correction chain
	function plotStatus(SC,Init,BPMData,CMData)

		REFPTS=1:length(Init.SC.RING);
		nREF = length(REFPTS);

		% DISTURBED LATTICE % % % % % % % % % % % % % % % % % % % % % % % %
		[ld,REAL.Tune,REAL.chromaticity] = atlinopt(Init.SC.RING,1e-3,REFPTS);
		REAL.beta                        = reshape([ld.beta],[2,nREF]);
		REAL.dispersion                  = reshape([ld.Dispersion],4,nREF);
		REAL.dispersion                  = REAL.dispersion([1 3],:);

		% IDEAL LATTICE % % % % % % % % % % % % % % % % % % % % % % % % % %
		% [ld,Ideal.Tune,Ideal.chromaticity] = atlinopt(IDEAL.Lattice,1e-3,REFPTS);
		[ld,Ideal.Tune,Ideal.chromaticity] = atlinopt(SC.IDEALRING,1e-3,REFPTS);
		Ideal.beta                         = reshape([ld.beta],[2,nREF]);
		Ideal.dispersion                   = reshape([ld.Dispersion],4,nREF);
		Ideal.dispersion                   = Ideal.dispersion([1 3],:);

		% CORRECTED LATTICE % % % % % % % % % % % % % % % % % % % % % % % %
		[ld,CorrModel.Tune,CorrModel.chromaticity] = atlinopt(SC.RING,1e-3,REFPTS);
		CorrModel.beta                             = reshape([ld.beta],[2,nREF]);
		CorrModel.dispersion                       = reshape([ld.Dispersion],4,nREF);
		CorrModel.dispersion                       = CorrModel.dispersion([1 3],:);

		BPMcal=[];CMcal=[];CMroll=[];fields={'H','V'};
		for nDim=1:2
			i=1;
			for ord=CMData.([fields{nDim} 'CMIndex'])'
				if nDim==1
					CMcal{nDim}(i)  = 1 + Init.SC.RING{ord}.CalErrorB(1);
				else
					CMcal{nDim}(i)  = 1 + Init.SC.RING{ord}.CalErrorA(1);
				end
				CMroll{nDim}(i) = Init.SC.RING{ord}.MagnetRoll(1) + Init.SC.RING{ord}.SupportRoll(1);
				i=i+1;
			end			
			BPMcal{nDim} = atgetfieldvalues(Init.SC.RING(BPMData.BPMIndex(BPMData.([fields{nDim} 'BPMIndex']))),'CalError',{1,nDim})';

		end

				
		% Plot BPM and Correctors
		fitCalCM   = {CMData.HCMKicks(:)'./1000./CMData.CMsteps{1}(:)'/2, CMData.VCMKicks(:)'./1000./CMData.CMsteps{2}(:)'/2};
		fitRollCM  = {-CMData.HCMCoupling',CMData.VCMCoupling'};
		fitCalBPM  = [BPMData.HBPMGain BPMData.VBPMGain]';
		fitRollBPM = [BPMData.VBPMCoupling-BPMData.HBPMCoupling]'/2;

		sPos = findspos(SC.RING,1:length(SC.RING))';

		titleBstr  = {'Hor. Beta','Ver. Beta'};titleBBstr = {'Hor. Dispersion','Ver. Dispersion'};
		figure(46),clf

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Beta beat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		for nDim=1:2
			bBeatInit(nDim,:) = 100*abs((Ideal.beta(nDim,:)-REAL.beta(     nDim,:))./Ideal.beta(nDim,:));
			bBeatCorr(nDim,:) = 100*abs((Ideal.beta(nDim,:)-CorrModel.beta(nDim,:))./Ideal.beta(nDim,:));

			subplot(2,6,6*(nDim-1)+1)
			plot(sPos,bBeatInit(nDim,:),'LineWidth',2);hold on
			plot(sPos,bBeatCorr(nDim,:),'LineWidth',2);
			xlabel('$S$ [m]');ylabel('$\Delta\beta/\beta_0$ [$\%$]');title(titleBstr{nDim})
			legend({sprintf('Init: %.1f $\\%%$',1E2*sqrt(mean(((Ideal.beta(nDim,:)-REAL.beta(nDim,:))./Ideal.beta(nDim,:)).^2))),...
				sprintf('Corr: %.1f $\\%%$',1E2*sqrt(mean(((Ideal.beta(nDim,:)-CorrModel.beta(nDim,:))./Ideal.beta(nDim,:)).^2)))},'Interpreter','Latex')
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Dispersion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		for nDim=1:2
			dBeatInit(nDim,:) = (Ideal.dispersion(nDim,:)-REAL.dispersion(nDim,:));
			dBeatCorr(nDim,:) = (Ideal.dispersion(nDim,:)-CorrModel.dispersion(nDim,:));

			subplot(2,6,6*(nDim-1)+2)
			plot(sPos,1E3*dBeatInit(nDim,:),'LineWidth',2);hold on
			plot(sPos,1E3*dBeatCorr(nDim,:),'LineWidth',2);
			xlabel('$S$ [m]');ylabel('$\Delta\eta$ [mm]');title(titleBBstr{nDim})
			legend({sprintf('Init: %.1f mm',1E3*sqrt(mean(dBeatInit(nDim,:).^2))),...
				sprintf('Corr: %.1f mm',1E3*sqrt(mean(dBeatCorr(nDim,:).^2)))})
		end

		for nDim=1:2

			subplot(2,6,3+(nDim-1)*6)
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% CM Calibration  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			plot(1E2*abs(1-CMcal{nDim}));hold on
			plot(1E2*abs(fitCalCM{nDim}-CMcal{nDim}))
			legend({sprintf('Init: $\\delta_{rms}$= %.1f$\\%%$',sqrt(mean((1E2*abs(1-CMcal{nDim})).^2))),sprintf('Corr: $\\delta_{rms}$= %.1f$\\%%$',sqrt(mean((1E2*abs(fitCalCM{nDim}-CMcal{nDim})).^2)))});
			ylabel('CM calibration error [$\%$]');xlabel('Index of CM');title('CM calibration error')

			subplot(2,6,4+(nDim-1)*6)
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% BPM Calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			x0 = 1E2*abs(BPMcal{nDim});
			x  = 1E2*abs(fitCalBPM(nDim,:)-(1+BPMcal{nDim}));
			plot(x0);hold on;plot(x)
			legend({sprintf('Init: $\\delta_{rms}$= %.1f$\\%%$',sqrt(mean(x0.^2))),sprintf('Corr: $\\delta_{rms}$= %.1f$\\%%$',sqrt(mean(x.^2)))});
			ylabel('BPM calibration error [$\%$]');xlabel('Index of BPM');title('BPM calibration error')

			subplot(2,6,5+(nDim-1)*6)
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% CM roll %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			x0 = 1E6*abs(CMroll{nDim});
			x  = 1E6*abs(fitRollCM{nDim}-CMroll{nDim});
			plot(x0);hold on;plot(x)
			legend({sprintf('Init: $\\delta_{rms}$= %.1furad',sqrt(mean(x0.^2))),sprintf('Corr: $\\delta_{rms}$= %.1furad',sqrt(mean(x.^2)))});
			ylabel('CM roll error [urad]');xlabel('Index of CM');title('CM roll error')

		end

		subplot(2,6,6)
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% BPM roll %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		x0 = 1E6*abs(atgetfieldvalues(Init.SC.RING(SC.ORD.BPM),'Roll')');
		x  = 1E6*abs(fitRollBPM-atgetfieldvalues(Init.SC.RING(SC.ORD.BPM),'Roll')');
		plot(x0);hold on;plot(x)
		legend({sprintf('Init: $\\delta_{rms}$= %.1furad',sqrt(mean(x0.^2))),sprintf('Corr: $\\delta_{rms}$= %.1furad',sqrt(mean(x.^2)))});
		ylabel('BPM roll error [urad]');xlabel('Index of BPM');title('BPM roll error')

		set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
		set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
		set(findall(gcf,'-property','FontSize'),'FontSize',14);
		set(gcf,'color','w');
		drawnow

	end
end
