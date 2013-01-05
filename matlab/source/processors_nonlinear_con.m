function [con,coneq] = processors_nonlinear_con(x,Rs,Ts,M,b,c,effp,par)
% Nonlinear constraints for optimization problem
% Including both order conditions and absolute monotonicity conditions
% 
% Used by pre_post_processors.m

%==========================================================================

%% Inequality constraints

% Extract arrays Rb,Tb from x
if par(1) == 1
    x1 = x(1:(Rs+1)*Rs/2);
    x2 = x(length(x1)+1:length(x1)+(Ts+1)*Ts/2);
else
    x1 = x(1:(Rs+3)*Rs/2);
    x2 = x(length(x1)+1:length(x1)+(Ts+3)*Ts/2);
end
[R,Rb] = unpack_RK(x1,Rs,par(1),0);
[T,Tb] = unpack_RK(x2,Ts,par(1),0);

Rr = -x(end-2); %Radius of absolute monotonicity of R
Tr = -x(end-1); %Radius of absolute monotonicity of T

%==========================================================================

% Absolute monotonicity conditions over Butcher coefficients

% Monotonicity of R
K1 = [R ; Rb'];
T1 = eye(size(R,1))+Rr*R;
con1 = K1/T1;
con2 = Rr*con1*ones(Rs,1)-ones(Rs+1,1);
con3 = par(2)*(sum(R,2)-1); % absicca less than one

% Monotonicity of T
K2 = [T ; Tb'];
T2 = eye(size(T,1))+Tr*T;
con4 = K2/T2;
con5 = Tr*con4*ones(Ts,1)-ones(Ts+1,1);
con6 = par(2)*(sum(T,2)-1); % absicca less than one

% Slack variable's inequalities: slack <= Rr and slack <= Tr 
con7 = x(end-2) - x(end);
con8 = x(end-1) - x(end);

con = [-con1(:); con2(:); con3(:); -con4(:); con5(:); con6(:); ...
    con7(:); con8(:)];

%==========================================================================

%% Equality constraints

% Order conditions
coneq = oc_processors(x,Rs,Ts,M,b,c,effp,par(1));

end 