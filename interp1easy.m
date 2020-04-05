function [YI XI] = interp1easy(Y,interpfac)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2013-11-13
% LAST MODIFIED: 2013-11-13
% Interpolate in one dimension using Matlab interp1 function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = 1:length(Y);
dxi = 1/round((length(Y))*interpfac-1);
XI = 0:dxi:1;
XI = XI*(length(Y)-1)+1;

YI = interp1(X,Y,XI,'spline');
