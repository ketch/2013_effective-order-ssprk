function dydt = stanescu1998(t,y)
% ODE system of two equations [Stanescu 1998]
%
% Used by myfun.m

%==========================================================================

dydt = [ 1/y(1) - y(2)*exp(t^2)/t^2 - t;
         1/y(2) - exp(t^2) - 2*t*exp(-t^2) ];
end