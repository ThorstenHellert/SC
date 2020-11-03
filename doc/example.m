


% === Setup enviroment

% Clear workspace, initialize the global variable `plotFunctionFlag`.
clear all
global plotFunctionFlag

% Set path to the AT function `atpath` and let AT set it's paths. Also set the
% path to the MML LOCO implementation and to the main SC folder which contains
% e.g. the function `SCinit`.
addpath('~/at/atmat');
atpath()
addpath('~/MML/applications/loco/');
addpath('~/sc');

% === Define lattice file

% Define a simple FODO lattice and print the summary.
QF   = atquadrupole('QF',...
	0.5, 1.2,...
	'PassMethod','StrMPoleSymplectic4RadPass',...
	'Energy',2.5E9);
QD   = atquadrupole('QD',...
	0.5,-1.2,...
	'PassMethod','StrMPoleSymplectic4RadPass',...
	'Energy',2.5E9);
SF   = atsextupole('SF',...
	0.1, 6.0487,...
	'PassMethod','StrMPoleSymplectic4RadPass',...
	'Energy',2.5E9);
SD   = atsextupole('SD',...
	0.1,-9.5203,...
	'PassMethod','StrMPoleSymplectic4RadPass',...
	'Energy',2.5E9);
BEND = atsbend('BEND',...
	1,2*pi/40,...
	'PassMethod','BndMPoleSymplectic4Pass',...
	'Energy',2.5E9);
RFC  = atrfcavity('RFCav','Energy',2.5E9);
D2   = atdrift('Drift',0.25);
D3   = atdrift('Drift',0.2);
MARK = @(name) atmarker(name,'IdentityPass');

cell = [{D2};{MARK('SectionStart')};...
	{MARK('GirderStart')};{BEND};{D3};{SF};{D3};{MARK('GirderEnd')};...
	{MARK('GirderStart')};{MARK('BPM')};{QF};{D2};{D2};{BEND};{D3};{SD};...
	{D3};{QD};{D2};{MARK('BPM')};{MARK('GirderEnd')};{MARK('SectionEnd')}];

RING = [{RFC};repmat(cell,20,1)];
RING = atsetcavity(RING,20e5,1,50);

atsummary(RING);

% === Initialize toolbox

% Initialize the SC toolbox with the previously defined lattice cell structure.
SC = SCinit(RING);

% === Register lattice in SC

% In the following section all relevant elements and error sources are
% registered in SC.

% Identify all BPMs in lattice structure and register them including
% uncertainties of the calibration factor, offset, roll, turn-by-turn noise and
% stored beam noise.
ords = SCgetOrds(SC.RING,'BPM');
SC = SCregisterBPMs(SC,ords,...
	'CalError',5E-2 * [1 1],... % relative
	'Offset',500E-6 * [1 1],... % [m]
	'Noise',10E-6 * [1 1],...   % [m]
	'NoiseCO',1E-6 * [1 1],...  % [m]
	'Roll',1E-3);               % [rad]

% Identify the QFs in the lattice structure and register them as horizontal
% corrector magnets with a limit of 1mrad and include uncertainties of the CM
% calibration factor, quadrupole strength error, magnet offset and magnet roll.
ords = SCgetOrds(SC.RING,'QF');
SC = SCregisterMagnets(SC,ords,...
	'HCM',1E-3,...                    % [rad]
	'CalErrorB',[5E-2 1E-3],...       % relative
	'MagnetOffset',200E-6 * [1 1],... % x and y, [m]
	'MagnetRoll',200E-6);             % [rad]

% Identify the QDs in the lattice structure and register them as vertical
% corrector magnets with a limit of 1mrad and include uncertainties of the CM
% calibration factor , quadrupole strength error, magnet offset and magnet
% roll.
ords = SCgetOrds(SC.RING,'QD');
SC = SCregisterMagnets(SC,ords,...
	'VCM',1E-3,...                    % [rad]
	'CalErrorA',[5E-2 0],...          % relative
	'CalErrorB',[0 1E-3],...          % relative
	'MagnetOffset',200E-6 * [1 1],... % x and y, [m]
	'MagnetRoll',200E-6);             % [rad]

% Identify the BENDs in the lattice structure and register them with a relative
% bending angle error and magnet offset and magnet roll.
ords = SCgetOrds(SC.RING,'BEND');
SC = SCregisterMagnets(SC,ords,...
	'BendingAngle',1E-3,...           % relative
	'MagnetOffset',200E-6 * [1 1],... % x and y, [m]
	'MagnetRoll',200E-6);             % [rad]

% Identify the SF&SD in the lattice structure and register them as skew
% quadrupole corrector magnets with a K value limit of 0.1 and include
% uncertainties of the skew quad calibration factor, sextupole strength error,
% magnet offset and magnet roll.
ords = SCgetOrds(SC.RING,'SF|SD');
SC = SCregisterMagnets(SC,ords,...
	'SkewQuad',0.1,...                 % [1/m]
	'CalErrorA',[0 1E-3 0],...         % relative
	'CalErrorB',[0 0 1E-3],...         % relative
	'MagnetOffset',200E-6 * [1 1],...  % x and y, [m]
	'MagnetRoll',200E-6);              % [rad]

%Identify the cavity in the lattice structure and register it including
%uncertainties for the frequency [Hz], voltage [V] and phase offset [m]
ords = findcells(SC.RING,'Frequency');
SC = SCregisterCAVs(SC,ords,...
	'FrequencyOffset',5E3,... % [Hz]
	'VoltageOffset',5E3,...   % [V]
	'TimeLagOffset',0.5);     % [m]

% Identify girder start and end ordinates in lattice structure and register
% them including uncertainties for the offset in x and y [m] and roll [rad]
ords = [SCgetOrds(SC.RING,'GirderStart');SCgetOrds(SC.RING,'GirderEnd')];
SC = SCregisterSupport(SC,...
	'Girder',ords,...
	'Offset',100E-6 * [1 1],... % x and y, [m]
	'Roll',200E-6);             % [rad]

% Identify section start and end ordinates in lattice structure and register
% them including uncertainties for the offset in x and y [m]
ords = [SCgetOrds(SC.RING,'SectionStart');SCgetOrds(SC.RING,'SectionEnd')];
SC = SCregisterSupport(SC,...
	'Section',ords,...
	'Offset',100E-6 * [1 1]); % x and y, [m]

% As a last registration step we define the 6x6 beam sigma matrix, random
% shot-to-shot injection variation and the uncertainty of the systematic
% injection errors, both in 6D.  We furthermore define the relative
% circumference uncertainty and the percentage of partcles which can be lost
% while still getting a proper BPM reading.
SC.INJ.beamSize = diag([200E-6, 100E-6, 100E-6, 50E-6, 1E-3, 1E-4].^2);

SC.SIG.randomInjectionZ = [1E-4; 1E-5; 1E-4; 1E-5; 1E-4; 1E-4]; % [m; rad; m; rad; rel.; m]
SC.SIG.staticInjectionZ = [1E-3; 1E-4; 1E-3; 1E-4; 1E-3; 1E-3]; % [m; rad; m; rad; rel.; m]

SC.SIG.Circumference = 2E-4; % relative
SC.BPM.beamLostAt    = 0.6;  % relative

% === Define lattice apertures

% In this section a simple aperture model is defined. The aperture radius of
% all drift spaces is 13mm, while an circular aperture is implemented in all
% magnets with a radius of 10mm. In order to create a `pinhole`, the 50th
% magnet is randomly choosen to get a small eliptical aperture.
for ord=SCgetOrds(SC.RING,'Drift')
	SC.RING{ord}.EApertures = 13E-3 * [1 1]; % [m]
end

for ord=SCgetOrds(SC.RING,'QF|QD|BEND|SF|SD')
	SC.RING{ord}.EApertures = 10E-3 * [1 1]; % [m]
end

SC.RING{SC.ORD.Magnet(50)}.EApertures = [6E-3 3E-3]; % [m]

% === Check registration

% In this section the SC registration is checked for consistency and the
% lattice is plotted.
SCsanityCheck(SC);

SCplotLattice(SC,'nSectors',10);

% === Apply errors

% The next step is to generate and apply an error set based on the previolusly
% defined uncertainties. The misalignments are plotted.
SC = SCapplyErrors(SC);

SCplotSupport(SC);

% === Setup correction chain

% At this point the parameters of the correction chain may be defined. In this
% example, we switch off the cavity and the sextupole magnets. Furthermore the
% 1 and 2-turn model trajectory response matrices are calcualted and a Tikhonov
% regularization with a regularization parameter of `50` is used to calculate
% the pseudo inverse of both matrices..
SC.RING = SCcronoff(SC.RING,'cavityoff');

sextOrds = SCgetOrds(SC.RING,'SF|SD');
SC = SCsetMags2SetPoints(SC,sextOrds,2,3,0,...
	'method','abs');

RM1 = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'nTurns',1);
RM2 = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'nTurns',2);

Minv1 = SCgetPinv(RM1,'alpha',50);
Minv2 = SCgetPinv(RM2,'alpha',50);

% Next, we define the number of particles per bunch, shots for averaging the
% BPM reading and number of turns and ensure turn-by-turn tracking mode. The
% noise level `eps` defines a stopping criteria for the feedback. Finally, we
% switch on the global plot falg and plot uncorrected beam trajectory.
SC.INJ.nParticles = 1;
SC.INJ.nTurns     = 1;
SC.INJ.nShots     = 1;
SC.INJ.trackMode  = 'TBT';

eps   = 1E-4; % Noise level

plotFunctionFlag = 0;

SCgetBPMreading(SC);

% === Start correction chain

% Run first turn feedback and apply correction if no error occured.
[CUR,ERROR] = SCfeedbackFirstTurn(SC,Minv1,'verbose',1);
if ~ERROR; SC=CUR; else; return; end

% Switch in 2-turn mode and get full 2-turn transmission by correcting the
% first three BPMs of the second turn to the corresponding readings in the
% first turn.
SC.INJ.nTurns = 2;

[CUR,ERROR] = SCfeedbackStitch(SC,Minv2,...
	'nBPMs',3,...
	'maxsteps',20,...
	'verbose',1);
if ~ERROR; SC=CUR; else; return; end

% Run trajectory feedback on 2-turn readings. Then create a period 1
% orbit by matching the second turn BPM readings to the first turn.
[CUR,ERROR] = SCfeedbackRun(SC,Minv2,...
	'target',300E-6,...
	'maxsteps',30,...
	'eps',eps,...
	'verbose',1);
if ~ERROR; SC=CUR; else; return; end

[CUR,ERROR] = SCfeedbackBalance(SC,Minv2,...
	'maxsteps',32,...
	'eps',eps,...
	'verbose',1);
if ~ERROR; SC=CUR; else; return; end

% In the following loop the sextupole magnets are ramped up in 5 steps
% and feedback is applied after each step.
for S = linspace(0.1,1,5)

	SC = SCsetMags2SetPoints(SC,sextOrds,2,3,S,...
		'method','rel');

	[CUR,ERROR] = SCfeedbackBalance(SC,Minv2,...
		'maxsteps',10,...
		'eps',eps,...
		'verbose',1);

	if ~ERROR; SC=CUR; end
end

% Switch off plotting every beam, switch the cavity on and plot initial
% phase space.
plotFunctionFlag = 0;

SC.RING = SCcronoff(SC.RING,'cavityon');

SCplotPhaseSpace(SC,...
	'nParticles',10,...
	'nTurns',100);

% The following block performs an rf phase and frequency correction in
% a loop and applies the corresponding correction step if no error
% occured.
for nIter=1:2
	% Perform RF phase correction.
	[deltaPhi,ERROR] = SCsynchPhaseCorrection(SC,...
		'nTurns',5,...      % Number of turns
		'nSteps',25,...     % Number of phase steps
		'plotResults',1,... % Final results are plotted
		'verbose',1);       % Print results
	if ERROR; error('Phase correction crashed');end

	% Apply phase correction
	SC = SCsetCavs2SetPoints(SC,SC.ORD.Cavity,...
			'TimeLag',deltaPhi,...
			'add');

	% Perform RF frequency correction.
	[deltaF,ERROR] = SCsynchEnergyCorrection(SC,...
		'range',40E3*[-1 1],... % Frequency range [kHz]
		'nTurns',20,...         % Number of turns
		'nSteps',15,...         % Number of frequency steps
		'plotResults',1,...     % Final results are plotted
		'verbose',1);           % Print results

	% Apply frequency correction
	if ~ERROR; SC = SCsetCavs2SetPoints(SC,SC.ORD.Cavity,...
			'Frequency',deltaF,...
			'add');
	else; return; end
end

% Plot final phase space and check if beam capture is achieved.
SCplotPhaseSpace(SC,'nParticles',10,'nTurns',1000);

[maxTurns,lostCount,ERROR] = SCgetBeamTransmission(SC,...
	'nParticles',100,...
	'nTurns',10,...
	'verbose',true);
if ERROR;return;end

% Beam capture achieved, switch to orbit mode for tracking. Calculate the orbit
% response matrix and the dispersion. Assume a beam based alignment procedure
% reduces the BPM offsets to 50um rms w.r.t. their neighbouring QF/QD magnets.
SC.INJ.trackMode = 'ORB';

MCO = SCgetModelRM(SC,SC.ORD.BPM,SC.ORD.CM,'trackMode','ORB');
eta = SCgetModelDispersion(SC,SC.ORD.BPM,SC.ORD.Cavity);

quadOrds = repmat(SCgetOrds(SC.RING,'QF|QD'),2,1);
SC       = SCpseudoBBA(SC,SC.ORD.BPM,quadOrds,50E-6);

% Run orbit feedback in a loop with decreasing Tikhonov regularization
% parameter `alpha` until no further improvment is achieved. Dispersion [m/Hz]
% is included and scaled by a factor of 1E8 to get the same magnitude as the
% orbit response [m/rad]
% matrix.
for	alpha = 10:-1:1
	% Get pseudo inverse
	MinvCO = SCgetPinv([MCO 1E8*eta],'alpha',alpha);

	% Run feedback
	[CUR,ERROR] = SCfeedbackRun(SC,MinvCO,...
		'target',0,...
		'maxsteps',50,...
		'scaleDisp',1E8,...
		'verbose',1);
	if ERROR;break;end

	% Calculate intial and final rms BPM reading.
	B0rms = sqrt(mean(SCgetBPMreading(SC ).^2,2));
	Brms  = sqrt(mean(SCgetBPMreading(CUR).^2,2));

	% Break if orbit feedback did not result in a smaller rms BPM reading
	if mean(B0rms)<mean(Brms);break;end

	% Accept new machine
	SC = CUR;
end

% === Perform LOCO based linear optics correction.

% The SC-LOCO interface is performed via a set of functions which are centrally
% stored in a pseudo-library `SClocoLib`.
% The first step is to setup the LOCO model
% ('setupLOCOmodel') from `SC`. Optional input
% arguments are passed to `LocoFlags`.  In this example dispersion is included
% in the evaluation and the horizontal and vertical weights are set to 100. The
% next function ('getBPMCMstructure') sets up the BPM and CM data structures.
% Again, optional arguments are passed to the corresponding structure, which
% allows to fit the CM and BPM calibration errros. Next, the orbit response
% matrix and the dispersion is measured ('getMeasurment')  using CM steps of
% 0.1mrad and an rf step of 1kHz, respectively. Finally the LOCO fit parameter
% structure is setup via 'setupFitparameters'. We start with all `QF` and `QD`
% quadrupoles which are individually powered and a strength variation of 1E-3
% and 1E-4 is used to calculate the derivatives, respectively.
CMstep = 1E-4; % [rad]
RFstep = 1E3;  % [Hz]

[LOCOmodel,LOCOflags,Init] = SClocoLib('setupLOCOmodel',SC,...
	'Dispersion','Yes',...
	'HorizontalDispersionWeight',.1E2,...
	'VerticalDispersionWeight',.1E2);

[BPMData,CMData] =  SClocoLib('getBPMCMstructure',SC,CMstep,...
	{'BPM','FitGains','Yes'},...
	{'CM','FitKicks','Yes'});

LOCOmeasData =  SClocoLib('getMeasurement',SC,CMstep,RFstep,SC.ORD.BPM,SC.ORD.CM);

FitParameters = SClocoLib('setupFitparameters',SC,Init.SC.RING,LOCOmodel,RFstep,...
	{SCgetOrds(SC.RING,'QF'),'normal','individual',1E-3},... % {Ords, normal/skew, ind/fam, deltaK}
	{SCgetOrds(SC.RING,'QD'),'normal','individual',1E-4});   % {Ords, normal/skew, ind/fam, deltaK}

% Run LOCO fit procedure in a loop and apply the lattice correction after each
% LOCO step, followed by an orbit correction step. After three iterations,
% include coupling (off-diagonal response matrix blocks) and the skew
% quadrupole correctors as LOCO fitparameters.
for n=1:6
	[~, BPMData, CMData, FitParameters, LOCOflags, LOCOmodel] = loco(LOCOmeasData,  BPMData,  CMData,  FitParameters,  LOCOflags,  LOCOmodel);

	SC = SClocoLib('applyLatticeCorrection',SC,FitParameters);

	SC = SClocoLib('applyOrbitCorrection',SC);

	SClocoLib('plotStatus',SC,Init,BPMData,CMData);

	if n==3
		LOCOflags.Coupling = 'Yes';

		FitParameters = SClocoLib('setupFitparameters',SC,Init.SC.RING,LOCOmodel,RFstep,...
			{SCgetOrds(SC.RING,'QF'),'normal','individual',1E-3},...
			{SCgetOrds(SC.RING,'QD'),'normal','individual',1E-4},...
			{SC.ORD.SkewQuad,'skew','individual',1E-3});
	end
end
