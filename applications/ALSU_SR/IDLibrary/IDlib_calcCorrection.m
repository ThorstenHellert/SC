function [deltaK,finalFOM] = IDlib_calcCorrection(SC,IDs,varargin)
% IDlib_calcCorrection
% ====================
%
% NAME
% ----
% IDlib_calcCorrection - Calculates the global ID compensation quadrupole values
%
% SYNOPSIS
% --------
% `[deltaK,finalFOM] = IDlib_calcCorrection(SC,IDs [,options])`
%
% DESCRIPTION
% -----------
% This function calcualtes the quadrupole setpont change required to compensate the linear optics
% pertubation by IDs. The IDs in the input AT lattice are closed.
%
% INPUT
% -----
% `SC`::
% 	The SC base structure including not closed IDs and without errors.
% `IDs`::
% 	Cell string with ID names.
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'IDstrength'` (`ones(size(IDs))`)::
%	Scaling factor for IDs (only for series of SBENDs)
% `'IDmap'` (`''`)::
%	Name of kick map to be loaded.
% `'regVal'` (`1`)::
%	Regularization value.
% `'regMethod'` (`'alpha'`)::
%	Regularization method.
% `'muWeighting'` (`100`)::
%	Weighting factor for phase advance vs. beta function error.
% `'magnetWeighting'` (`[]`)::
%	Weighting factors for individual magnets.
% `'dKstop'` (`1E-3`)::
%	Quadrupole setpoint change stopping criteria.
% `'quadOrds'` (`SCgetOrds(SC.RING,'^QD1|^QF[1-2]|^QF[4-6]|BEND2|BEND3')`)::
%	Quadrupole ordinates for ID compensation.
% `'evalOrds'` (`SCgetOrds(SC.RING,'SF|SD')`)::
%	Ordinates for figure of merit calculation.
% `'verbose'` (`0`)::
%	If true debug information is printed.
% `'plot'` (`0`)::
%	If true results are plotted.
%
% RETURN VALUE
% ------------
% `deltaK`::
% 	Vector of quadrupole setpoint changes.
% `finalFOM`::
% 	Value of final figure of merrit.
%
% SEE ALSO
% --------
% *IDLib_closeID*, *IDlib_includeIDs*, *IDLib_getResponseMatrix*, *IDlib_applyCorrection*, 


	% Parse optional arguments
	p = inputParser;
	addOptional(p,'IDstrength',ones(size(IDs)));
	addOptional(p,'IDmap',[]);
	addOptional(p,'regVal',1);
	addOptional(p,'regMethod','alpha');
	addOptional(p,'muWeighting',100);
	addOptional(p,'magnetWeighting',[]);
	addOptional(p,'dKstop',1E-3);
	addOptional(p,'quadOrds',SCgetOrds(SC.RING,'^QD1|^QF[1-2]|^QF[4-6]|BEND2|BEND3'));
	addOptional(p,'evalOrds',SCgetOrds(SC.RING,'SF|SD'));
	addOptional(p,'verbose',0);
	addOptional(p,'plot',0);
	parse(p,varargin{:});
	par = p.Results;

	
	% Check if magnet weighting should not be applied
	if isempty(par.magnetWeighting)
		par.magnetWeighting = ones(size(par.evalOrds));
	end
	
	% Get reference values
	IdealMat = IDLib_getResponseMatrix(SC,par.quadOrds,par.evalOrds,'muWeighting',par.muWeighting,'magnetWeighting',par.magnetWeighting);
	RefVal   = IDLib_getResponseMatrix(SC,par.quadOrds,par.evalOrds,'refOnly',1,'muWeighting',par.muWeighting,'magnetWeighting',par.magnetWeighting);
	
	% Get quadrupole startpoints
	for nMagnet=1:length(par.quadOrds)
		startPoints(nMagnet) = SC.RING{par.quadOrds(nMagnet)}.SetPointB(2);
	end
	
	% For plotting
	[out,~,~] = atlinopt(SC.RING,0,1:length(SC.RING));
	beta0 = vertcat(out.beta);
	sPos = findspos(SC.RING,1:length(SC.RING))';
	
	% Close IDs
	for nID=1:length(IDs)
		SC = IDLib_closeID(SC,IDs{nID},'IDstrength',par.IDstrength(nID),'IDmap',par.IDmap{nID},'PassMethod','ThinEPU2Pass');
	end
	
	% For plotting
	[out,~,~] = atlinopt(SC.RING,0,1:length(SC.RING));
	beta1 = vertcat(out.beta);
	
	
	
	
	fprintf('Regularization %.1f\n',par.regVal)
	
	% Get pseudo inverse
	Minv     = SCgetPinv(IdealMat,par.regMethod,par.regVal,'plot',0,'damping',.1);
	
	% Number of feedback steps
	for nStep=1:50
		
		% Get current figure of merrit value
		RealVal = RefVal - IDLib_getResponseMatrix(SC,par.quadOrds,par.evalOrds,'refOnly',1,'muWeighting',par.muWeighting,'magnetWeighting',par.magnetWeighting);
		
		% Break if something went wrong
		if any(isnan(RealVal))
			break
		end
		
		% Calculate correction step
		dK = Minv * RealVal;
		
		if par.verbose;fprintf('Std(dK) = %.1e\n',std(dK));end
		
		% Apply correction step
		SC = SCsetMags2SetPoints(SC,par.quadOrds,2,2,-dK,'method','add','dipCompensation',1);
		
		% Plot current state
		if par.plot;plotCurrentStep();end
		
		% Check if stopping criteria is met
		if std(dK)<par.dKstop
			break
		end
		
	end
	
	% Get final FOM value
	finalFOM = std(RefVal - IDLib_getResponseMatrix(SC,par.quadOrds,par.evalOrds,'refOnly',1,'muWeighting',par.muWeighting,'magnetWeighting',par.magnetWeighting));
	
	% Get final correction step
	for nMagnet=1:length(par.quadOrds)
		deltaK(nMagnet) = SC.RING{par.quadOrds(nMagnet)}.SetPointB(2) - startPoints(nMagnet);
	end
	
	if par.verbose;fprintf('FOM = %.1e\n',finalFOM);end
	
	
	

	
	function plotCurrentStep
		
		[out,~,~] = atlinopt(SC.RING,0,1:length(SC.RING));
		beta2 = vertcat(out.beta);
		
		figure(42);clf
		for nDim=1:2
			subplot(4,1,nDim);hold on
			plot(sPos,beta1(:,nDim),'LineWidth',2)
			plot(sPos,beta2(:,nDim),'LineWidth',2)
			plot(sPos,beta0(:,nDim),'k','LineWidth',2)
			xlabel('s [m]');ylabel('Beta [m]');legend({'Dist','Corr','Ideal'});
			subplot(4,1,2+nDim);hold on
			plot(sPos,beta1(:,nDim)-beta0(:,nDim),'LineWidth',2)
			plot(sPos,beta2(:,nDim)-beta0(:,nDim),'LineWidth',2)
			xlabel('s [m]');ylabel('$\Delta$ Beta [m]');legend({'Dist','Corr'});
		end
		drawnow
	end
	

end
