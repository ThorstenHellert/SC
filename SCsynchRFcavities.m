function [deltaPhis,maxTurns,ERROR] = SCsynchRFcavities(SC,varargin)
% SCsynchRFcavities
% ======================
%
% NAME
% ----
% SCsynchRFcavities - Returns the phases that maximize the number of turns.
%
% SYNOPSIS
% --------
% `[deltaPhis,maxTurns,ERROR] = SCsynchRFcavities(SC [, options])`
%
%
% DESCRIPTION
% -----------
%
%  
% This function is useful when commissioning multiple cavities in TBT mode.
% It scans the phase interval [-pi,pi] stepwise one cavity at the time and
% returns the phases achieving the maximum number of turns.
%
%
% INPUTS
% ------
% `SC`::
%   The SC base structure
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'cavOrd'` (`SC.ORD.Cavity`)::
%	Ordinate of evaluated cavity
% `'nSteps'` (`20`)::
%	Number of phase steps to be evaluated
% `'nTurns'` (`SC.INJ.nTurns`)::
%	Number of turns to check beam transmission.
% `'verbose'` (`0`)::
%	If 1, status is printed per cavity.
%	If 2, status is printed per cavity and step.
% `'plot'` (`0`)::
%	If true, the max. number of turns vs the phase offset is plotted.
%
% RETURN VALUES
% -------------
% `deltaPhis`::
%   Row array of length equal to the number of cavities.
%   It contains the phase to be added to cavity field 'TimeLag'.
% `maxTurns`::
%   N by M matrix for N cavities and M steps.
%   It contains the maximum number of turns per cavity and setpoint.
% `ERROR`::
%   N by M matrix for N cavities and M steps.
%   It contains the ERROR flag from SCgetBeamTransmission.
%
% ERRORS
% ------
% `False`::
% 	Beam survives all turns. See SCgetBeamTransmission.
%
% EXAMPLE
% -------
% [Phs,~,~] = SynchRFCavities(SC,'nturns',nturns,'nSteps',40);
% SC = SCsetCavs2SetPoints(SC,SC.ORD.Cavity(1),'TimeLag',Phs(1),'add');
% SC = SCsetCavs2SetPoints(SC,SC.ORD.Cavity(2),'TimeLag',Phs(2),'add');
%
% SEE ALSO
% --------
% *SCsetCavs2SetPoints*, *SCgetBPMreading*, *SCsynchPhaseCorrection*,
% *SCgetBeamTransmission*

%  Author     Z. Marti  at ALBA CELLS, 2024nov
%  Edited by  O. Blanco at ALBA CELLS, 2025feb

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Input check

	% Parse input
    p = inputParser;
    addOptional(p,'cavOrd',SC.ORD.Cavity);
    addOptional(p,'nSteps',20);
    addOptional(p,'nturns',SC.INJ.nTurns);
    addOptional(p,'verbose',0);
    addOptional(p,'plot',0);
    parse(p,varargin{:});
    par = p.Results;

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Initialization
    
    % number of turns
    nturns = par.nturns;

    % plotting
    doplot = par.plot;

    % verbose mode
    doverbose = par.verbose;

    % get the cavities Ords
    cavOrds = par.cavOrd;

    % get the total number of cavities
    ncav = numel(cavOrds);

    % RF wavelength [m]
    meanfreq = mean(atgetfieldvalues(SC.IDEALRING,cavOrds,'Frequency'));
    lambda = 299792458/meanfreq;

    % number of phase setpoints
    nSteps = par.nSteps;

    % Define phase scan range
    lambdaTestVec = 1/2 * lambda * linspace(-1,1,nSteps);

    % store the maximum number of turn per cavity and phase setpoint
    maxTurns  = zeros(ncav,nSteps);

    % store the lost count from the cavities and phase setpoints
    lost = cell(ncav,nSteps);

    % store the transmission errors from 
    ERROR = zeros(nSteps,nSteps);

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Main script
    
    % count the number of cavities
    if doverbose
        fprintf('Synchronization of cavities ... \n');
    end

    % loop over the cavities
    for ii=1:ncav

        % print cavity index and order
        if doverbose == 1
            fprintf("Cavity %d of %d, ord %d\n",ii,ncav,cavOrds(ii));
        end

        % turn off all cavities
        for kk=1:ncav
            SC.RING{cavOrds(kk)}.PassMethod = 'IdentityPass';
        end

        % switch on one cavity
        SC.RING{cavOrds(ii)}.PassMethod = 'RFCavityPass';

        % loop nSteps times over the phase and check transmission turns
        for jj=1:nSteps

            % print cavity index and order
            if doverbose == 2
                fprintf("Cavity %d of %d, ord %d, step %d of %d\n", ...
                        ii,ncav,cavOrds(ii),jj,nSteps);
            end

            % offset the initial the cavity phase
            tmpSC = SCsetCavs2SetPoints(SC,cavOrds(ii), ...
                                           'TimeLag',lambdaTestVec(jj),...
                                           'add');

            % check the maximum number of turns achieved by the beam
            [ maxTurns(ii,jj), ...
              lost{ii,jj}, ...
              ERROR(ii,jj) ] = SCgetBeamTransmission(tmpSC, ...
                                                     'nturns',nturns);
        end

    end

    % get the setpoint that achieves maximum transmission per cavity
    [~, mi] = max(maxTurns,[],2);

    % get the best phase from the setpoint vector
    deltaPhis = lambdaTestVec(mi);

    % plot
    if doplot
        figure(89);
        clf(89);
        cmp = hsv(ncav);
        hold on;
        for idxcav = 1:ncav
            plot(lambdaTestVec,maxTurns(idxcav,:), ...
                '-o', ...
                'DisplayName',sprintf('Cavity %d',idxcav), ...
                'Color',cmp(idxcav,:));
            plot(lambdaTestVec(mi(idxcav)),maxTurns(idxcav,mi(idxcav)), ...
                'x', ...
                'Color',cmp(idxcav,:), 'MarkerSize', 24, ...
                'HandleVisibility','off');
        end
        xlabel('$\Delta \phi$ [m]');
        ylabel('Max. Number of Turns');
        legend;
        set(findall(gcf,'-property','FontSize'),'FontSize',18);
	    set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
	    set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
	    set(gcf,'color','w');
    end

    % plot 2D, valid only for two cavities
    % only useful when scanning one vs another
    % if doplot 
    %     if ncav ~= 2
    %         fprintf("Plotting is only possible for 2 cavities.\n");
    %     else
    %         figure(89);
    %         pcolor( ...
    %             lambdaTestVec(:) * ones(1,nSteps), ...
    %             ones(nSteps,1) * lambdaTestVec, ...
    %             maxTurns);
    %         shading flat;
    %         colorbar;
    %         title('Max. turns');
    %         hold all;
    %         plot(lambdaTestVec, lambdaTestVec, '-r');
    %         xlabel('CT lag CAV 1 [m]');
    %         ylabel('CT lag CAV 2 [m]')
    %     end
    % end

    if doverbose == 1 || doverbose == 2
        fprintf(" ... done.\n")
    end

end