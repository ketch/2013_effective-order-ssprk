function [alpha,beta,A,b,c,r,ceff] = load_method(s,method)
% Returns the optimsl SSP Runge-Kutta method and SSP coefficient for 
% particular cases. 
% For these cases, formulas of the method's coefficients are known.
% Cases: 1. ESSPRK(v,2,2), v > 1
%        2. ESSPRK(v^2,3,3), v > 2
%        3. ESSPRK(v^2+1,4,2), v > 2
%
% Used by effective_ssp.m

%==========================================================================

if strcmp(method,'ESSPRK2')
    r = s-1;
    
    % alpha matrix
    alpha = diag(ones(s,1),-1);
    alpha(end,s) = (s-1)/s;
    alpha(end,1) = 1/s;
    alpha = alpha(:,1:end-1); % fixing dimensions
    
    % beta matrix 
    beta = diag(1/(s-1)*ones(s,1),-1);
    beta(end,s) = 1/s;
    beta = beta(:,1:end-1); % fixing dimensions
    
    % Butcher tableau and effective SSP coefficient
    [A,b,c] = shuosher2butcher(alpha,beta);
    ceff = r/s;
    
elseif strcmp(method,'ESSPRK3')
    v = sqrt(s);
    r = v^2-v;
    
    % alpha matrix
    alpha = diag(ones(s,1),-1);
    alpha(v*(v+1)/2+1,v*(v+1)/2) = (v-1)/(2*v-1);
    alpha(v*(v+1)/2+1,(v-1)*(v-2)/2+1) = v/(2*v-1);
    alpha = alpha(:,1:end-1); % fixing dimensions
    
    % beta matrix 
    beta = alpha/r;
    beta(v*(v+1)/2+1,(v-1)*(v-2)/2+1) = 0;
    
    % Butcher tableau and effective SSP coefficient
    [A,b,c] = shuosher2butcher(alpha,beta);
    ceff = r/s;
    
elseif strcmp(method,'ESSPRK4')
    v = sqrt(s-1);
    r = v^2-v;
    
    % alpha matrix
    alpha = zeros(s+1,s+1);
    
    alpha(v^2-2*v+4,(v-2)^2) = (v^2-1 - sqrt(v^3-3*v^2+v+1))/(4*v^2-6*v+2); 
    %also (v^2-1 + sqrt(v^3-3*v^2+v+1))/(4*v^2-6*v+2)
    alpha(end,1) = 2/(s*((v-1)^2+1));
    alpha(end,v^2-2*v+2) = 1 - alpha(end,1) - ...
        (r*(v-1))/(s*(2*v-1)*(1-alpha(v^2-2*v+4,(v-2)^2)));
    
    % Construct alpha matrix
    bigdiag = 1 - sum(alpha,2);
    alpha = alpha + diag(bigdiag(2:end),-1);
    alpha = alpha(:,1:end-1); % fixing dimensions
    
    % beta matrix 
    beta = alpha/r;
    beta(end,1) = 0;
    
    % Butcher tableau and effective SSP coefficient
    [A,b,c] = shuosher2butcher(alpha,beta);
    ceff = r/s;
    
end

end
    