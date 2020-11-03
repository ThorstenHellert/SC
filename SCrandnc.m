function out=SCrandnc(c,varargin)
% SCrandnc
% ========
%
% NAME
% ----
% SCrandnc - Truncated Gaussian random number generator
%
% SYNOPSIS
% --------
% `out = SCrandnc(c, [, options])`
%
%
% DESCRIPTION
% -----------
% Random number generator (Gaussian distribution with mean 0 and sigma of 1). This function calls
% the built in Matlab function randn but allows to specify a cut parameter `c`. Only numbers between
% `-c` and `+c` will be generated.
%
% INPUTS
% ------
% `c`:: cutoff value
%
%
% OPTIONS
% -------
% `size`:: Either [`M`, `N`] array defining the size of the output array or two values `M`, `N`. See examples.
%
% RETURN VALUES
% -------------
% `output`::   [`M` x `N`] array of random numbers with a Gaussian distribution
%
%
% EXAMPLES
% --------
% Generates a random number `out` between [-2,2] and with sigma=1.
% ------------------------------------------------------------------
% out = SCrandnc(2);
% ------------------------------------------------------------------
%
% Generates random numbers `out` between [-2,2] and with sigma=1 with size(out)=size(A).
% ------------------------------------------------------------------
% out = SCrandnc(2,size(A));
% ------------------------------------------------------------------
%
% Generates random numbers `out` between [-2,2] and with sigma=1 with size(out)=[1,10].
% ------------------------------------------------------------------
% out = SCrandnc(2,1,10);
% ------------------------------------------------------------------


% Check for normalization flag
if any(strcmp(varargin,'normalize'))
	normalize = 1;
	varargin(strcmp(varargin,'normalize')) = [];
else
	normalize = 0;
end

% Define size of output
if isempty(varargin)
	m = 1;
	n = 1;
elseif length(varargin)==1
	m = varargin{1}(1);
	n = varargin{1}(2);
elseif length(varargin)==2
	m = varargin{1};
	n = varargin{2};
else
	error('Too much input arguments.')
end

% Generate random numbers
out = randn(m,n);

% Find outlier index
outindex = find(abs(out)>abs(c));

% Fill outliers with new values
while ~isempty(outindex)
    out(outindex) = randn(size(outindex));
    outindex = find(abs(out)>abs(c));
end

% Normalize output to sigma=1
if normalize
	warning('Not yet implemented.')
end
