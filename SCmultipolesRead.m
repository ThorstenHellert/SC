function [AB, order, type] = SCmultipolesRead(fname)
% SCmultipolesRead
% ================
%
% NAME
% ----
% SCmultipolesRead - Reads multipole tables
%
% SYNOPSIS
% --------
% `[AB, order, type] = SCmultipolesRead(fname)`
%
%
% DESCRIPTION
% -----------
% This functions reads multipole coefficients from the text-file `fname`, see
% below for an example.
% This file needs to contain all orders up to the maximum order; no gaps!
% Lines starting with `#` are ignored.
% Fields are delimited by tabs.
%
%
% INPUT
% -----
% `fname`::
%	File to read.
%
% RETURN VALUES
% -------------
% `AB`::
%	`[2,maxorder]`-array containing the multipole coefficients
% `order`, `type`::
%	If the multipole table contains an exact `1.0`, it is assumed
%	this table is scaled to that coefficient. In that case
%	`order` gives the `n` of that coefficent and `type` determines
%	whether it is the normal (2) or skew (1) component.
%
% EXAMPLE FILE
% ------------
% ---------------------------------
% # Taken from v20r_SysMultipoles_R1.xlsx (table DIPA, AT conv)
% # This is a comment.
% # n	PolynomA(n)		PolynomB(n)
% 0	+0.000000000E+00	-1.712194700E-01
% 1	+0.000000000E+00	+2.997864116E+00
% 2	+0.000000000E+00	+1.121256852E+00
% 3	+0.000000000E+00	-5.000867669E+01
% 4	+0.000000000E+00	-4.233346976E+03
% 5	+0.000000000E+00	-6.454146352E+05
% 6	+0.000000000E+00	-4.014924846E+07
% 7	+0.000000000E+00	-8.862779057E+09
% 8	+0.000000000E+00	+2.614134332E+12
% 9	+0.000000000E+00	-5.382882472E+14
% ---------------------------------


f=fopen(fname,'r');
	tab=cell2mat(textscan(f,'%f%f%f','Delimiter','\t','CommentStyle','#'));
	fclose(f);
	if(size(tab,2)~=3) error('Incorrect table size.'); end;
	AB = tab(:,2:end);
	idx=find(AB==1);
	if nargout>1 && length(idx)~=1;
		warning('Nominal order could not be (uniquely) determined. Continuing with idx=1.');
		idx=1;
	end;
	[order,type] = ind2sub(size(AB),idx);
	if type>2; error('Ill-defined magnet type.'); end;
end
