function [dbz] = dbzero(mat,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% CREATED: 2012-01-23
% LAST MODIFIED: 2020-09-02
% db scale with max at zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  optargin = size(varargin,2);
  if(optargin==1)
    mx=db(varargin{1});
  else
    mx=maxmax(db(mat));
  end
dbz=db(mat);
dbz=dbz-mx;
