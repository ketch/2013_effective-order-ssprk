function [con,coneq] = ssp_nonlinear_con(x,s,p,effp,par)
% Nonlinear constraints for optimization problem
% Including both order conditions and absolute monotonicity conditions
% 
% Used by effective_ssp.m

%==========================================================================

%% Inequality constraints

% Radius of absolute monotonicity
r = -x(end);

if par(2) == 0
    % Extract coefficient arrays from x 
    [A,b] = unpack_RK(x,s,par(1),par(2));
    
    % Absolute monotonicity conditions over Butcher coefficients
    K = [A; b'];
    T = eye(s)+r*A;
    
    con1 = K/T;
    con2 = r*con1*ones(s,1) - ones(s+1,1);
else
    % Extract coefficient arrays from x 
    [alpha,beta] = unpack_RK(x,s,par(1),par(2));
    
    % Absolute monotonicity conditions over Shu-Osher coefficients
    con1 = x(1:end-1);
    con2 = r*beta*ones(s,1) - ones(s+1,1);
    
    A = shuosher2butcher(alpha,beta);
end

con3 = par(2)*(sum(A,2) - ones(s,1)); % absicca less than one
con = [-con1(:); con2(:); con3(:)];

%==========================================================================

%% Equality constraints

% Order conditions
coneq = oc_effective_ssp(x,s,p,effp,par(1),par(2));

end
