function dydt = krogh(t,y)
% Euler equations of a rigid body without external forces
%
% Used by myfun.m

%==========================================================================

dydt = [ y(2)*y(3);
         -y(1)*y(3);
         -0.51*y(1)*y(2) ];

end

