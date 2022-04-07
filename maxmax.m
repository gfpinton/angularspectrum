function [mat] = maxmax(mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2012-01-23
% LAST MODIFIED: 2013-11-13
% max of a matrix with arbitrary dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nd = length(size(mat));

for n=1:nd
  mat = max(mat);
end
