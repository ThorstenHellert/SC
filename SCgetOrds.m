function ords = SCgetOrds(RING,rx,varargin)
% SCgetOrds
% =========
%
% NAME
% ----
% SCgetOrds - get ordinates of elements from FamName
%
% SYNOPSIS
% --------
% `ords = SCgetOrds(RING, rx [, options])`
%
%
% DESCRIPTION
% -----------
% SCgetOrds produces the ordinates (indices) of elements in `RING` whose
% `FamName` matches the regular expression(s) `rx`.  `rx` can be a single regex
% or a cell-array containing multiple regexes.
%
%
% RETURN VALUE
% ------------
% `ords`::
%	If `rx` is a single string `ords` is a single array; if `rx` is a
%	cell-array `ords` is a cell-array containing the respective ordinates.
%
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'verbose'` (0)::
%	If true additional information is printed.
%
% EXAMPLES
% --------
% Find all lattice elements in `SC.RING` beginning with `QF` and `QD` and return 
% them in a single array `ords`.
% ------------------------------------------------------------------
% ords = SCgetOrds(SC.RING,'^QF|^QD');
% ------------------------------------------------------------------
%
% Find all lattice elements in `SC.RING` named `QF` and `QD` and return 
% each of them in a cell array `ords`.
% ------------------------------------------------------------------
% ords = SCgetOrds(SC.RING,{'^QF$','^QD$'});
% ------------------------------------------------------------------


	p = inputParser;
	addOptional(p,'verbose',0);
	parse(p,varargin{:});
	verbose = p.Results.verbose;

	ords = [];

	if iscellstr(rx)
		for r=rx
			ords{end+1} = SCgetOrds(RING,r{1},varargin{:});
		end
		return
	end

	for ord=1:length(RING)
		if regexp(RING{ord}.FamName,rx)
			ords(end+1)=ord;
			if verbose; fprintf('Matched: %s\n',RING{ord}.FamName);
		end
	end
end
