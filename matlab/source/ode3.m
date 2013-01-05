function dydt = ode3(t,y)
% ODE system of 3 equations
%
% Used by myfun.m

%==========================================================================

dydt = [ -y(1);
         y(1) - y(2)^2;
         y(2)^2 ];
end