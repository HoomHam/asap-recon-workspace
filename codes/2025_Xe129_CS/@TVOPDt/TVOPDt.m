function  res = TVOPDt()

%res = TVOP()
%
% Implements a spatial finite-differencing operator.
%
% (c) Michael Lustig 2007

res.adjoint = 0;
res = class(res,'TVOPDt');

