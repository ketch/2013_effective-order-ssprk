function dydt = ode2(t,y)
% ODE system of 2 equations
%
% Used by myfun.m

%==========================================================================

dydt = [ y(1);
         2*y(2) ];
end