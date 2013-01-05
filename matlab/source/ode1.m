function dydt = ode1(t,y)
% ODE system of 2 equations
%
% Used by myfun.m

%==========================================================================

dydt = [ 2*(y(1)-y(1)*y(2));
         -(y(2)-y(1)*y(2)) ];
end