function dydt = vdp(t,y)
% van der Pol's equation
%
% Used by myfun.m

%==========================================================================

dydt = [ y(2);
         2.*(1-y(1)^2)*y(2)-y(1) ];
     
end
