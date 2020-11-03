function Minv = SCgetPinv(M,varargin)
% SCgetPinv
% =========
%
% NAME
% ----
% SCgetPinv - Calculates the pseudo inverse of a matrix
%
% SYNOPSIS
% --------
% `Minv = SCgetPinv(M [, options])`
%
%
% DESCRIPTION
% -----------
% Calculats the pseudo-inverse `Minv` of the input matrix `M` based on a singular value
% decomposition. Optional parameters define the number of singular values to be set to zero, a
% global scaling factor applied to all singular values or the Tikhonov regularization parameter.
%
% INPUTS
% ------
% `M`::  Input matrix
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'N'` (0)::
%   Number of singular values to be set to zero (starting from smallest value)
% `'damping'` ([])::
%   If not empty, all singular values are scaled by this factor
% `'alpha'` ([])::
%   If not empty, Tikhonov regularization is applied with regularization parameter `alpha`
% `'plot'` (0)::
%   If true, initial and inverse singular value spectrum is plotted
%
% RETURN VALUE
% ------------
% `Minv`::
% 	Pseudo-inverse of the input matrix.
%
%
% EXAMPLES
% --------
% Calculates the svd-based pseudo-inverse
% ------------------------------------------------------------------
% Minv = SCgetPinv(M);
% ------------------------------------------------------------------
%
% Calculates the svd-based pseudo-inverse while cutting off the last 10 singular values.
% ------------------------------------------------------------------
% Minv = SCgetPinv(M,'N',10);
% ------------------------------------------------------------------
%
% Calculates the svd-based pseudo-inverse while cutting off the last 10 singular values and using a
% Tikhonov regularization with regularization parameter of 10.
% ------------------------------------------------------------------
% Minv = SCgetPinv(M,'N',10,'alpha',1);
% ------------------------------------------------------------------
%
% SEE ALSO
% --------
% *SCgetModelRM*


	% Parse input
	p = inputParser;
	addOptional(p,'N',0);
	addOptional(p,'alpha',[]);
	addOptional(p,'damping',[]);
	addOptional(p,'plot',0);
	parse(p,varargin{:});
	par = p.Results;

	% Singular value decomposition
	[U,S,V] = svd(M);

	% Get singular values
	SVs = diag(S);


	D = zeros(size(S));
	% Perform Tikhonov regularization
	if ~isempty(par.alpha)
		D(1:length(SVs),1:length(SVs)) = diag(SVs ./ (SVs .* SVs + par.alpha^2));
	else
		D(1:length(SVs),1:length(SVs)) = diag(1./SVs);
	end

	% Apply damping factor
	if ~isempty(par.damping)
		D = par.damping * D;
	end

	% Cut singular values
	if par.N~=0
		keep = length(SVs)-par.N;
		D(keep+1:end,keep+1:end) = 0;
	end

	% Calculate final pseudo inverse.
	Minv = V * D' * U';

	if par.plot
		figure(66);
		subplot(1,2,1);
		semilogy(diag(S)/max(diag(S)),'o--');
		xlabel('Number of SV');ylabel('$\sigma/\sigma_0$');
		subplot(1,2,2);
		plot(diag(S).*diag(D),'o--')
		xlabel('Number of SV');ylabel('$\sigma * \sigma^+$');set(gca,'yaxislocation','right')
		set(findall(gcf,'-property','FontSize'),'FontSize',18);
		set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
		set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
		set(gcf,'color','w');
		drawnow
	end
end



