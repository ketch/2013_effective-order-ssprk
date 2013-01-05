function dydt = nonlinear7(t,y)
% Single ODE equations
%
% Used by myfun.m

%==========================================================================

dydt = (y - t)/(y + t);

end