function [array1,array2] = unpack_RK(x,s,explicit,shuoshercoefficients)
% Extracts the coefficient arrays from the optimization vector
%
% Used by effective_ssp.m and oc_effective_ssp.m

%==========================================================================

if shuoshercoefficients == 0
    %% Unpack optimization vector to Butcher matrix A and vector b
    
    % Construct matrix A and vector b
    A = zeros(s);
    if explicit == 1
        for i=2:s
            A(i,1:i-1) = x((i-2)*(i-1)/2+1:i*(i-1)/2);
        end
        b = x(s*(s-1)/2+1:s*(s-1)/2+s)';
    else 
        for i=1:s
            A(i,1:i) = x(i*(i-1)/2+1:i*(i+1)/2);
        end
        b = x(s*(s+1)/2+1:s*(s+1)/2+s)';
    end
    
    % Outputs
    array1 = A;
    array2 = b;
else
    %% Unpack optimization vector to Shu-Osher matrices alpha and beta
    
    % Construct beta matrix
    beta = zeros(s+1,s);
    if explicit == 1
        for i=2:s+1
            beta(i,1:i-1) = x((i-2)*(i-1)/2+1:i*(i-1)/2);
        end
    else
        for i=1:s
            beta(i,1:i) = x(i*(i-1)/2+1:i*(i+1)/2);
        end
        beta(s+1,1:s) = x(s*(s+1)/2+1:s*(s+1)/2+s);
    end
    
    % Construct alpha matrix
    r = -x(end);
    alpha = r*beta;
    
    % Consistency for alpha matrix in the Shu-Osher form
    for i=2:s+1
        alpha(i,1)=1-sum(alpha(i,2:s));
    end
    
    % Outputs
    array1 = alpha;
    array2 = beta;
end

end

% % ESSP(s^2+1,4,2)
% if strcmp(essp_method,'special')
%     v = sqrt(s-1);
%     alpha = zeros(s+1,s+1);
%     
%     % plugging optimization variables in alpha matrix
%     alpha(v^2-2*v+4,(v-2)^2) = x(1);
%     alpha(end-3,v^2-2*v-1) = 0;%x(2);
%     alpha(end-1,v^2-2*v+1) = x(3);
%     alpha(end,v^2-2*v+2) = x(4);
%     alpha(end,1) = 2/(s*((v-1)^2+1));
%     
%     % Construction alpha matrix
%     big_diag = 1 - sum(alpha,2);
%     alpha = alpha + diag(big_diag(2:end),-1);
%     alpha = alpha(:,1:end-1); % fixing dimensions
%     
%     % Construct beta matrix
%     r = -x(end); 
%     alpha(end,1) = 0;
%     beta = alpha/r;
%     
%     % Consistency for alpha matrix in the Shu-Osher form
%     for i=2:s+1
%        alpha(i,1)=1-sum(alpha(i,2:s));
%     end
%     
% elseif strcmp(essp_method,'special2')
%     v = sqrt(s-1);
%     %r = v^2-v;
%     r = -x(end); 
%     %x(1) = 0; x(3) = 0;                    % 10 stages, v=3
%     %x(1:4) = 0; x(6:7) = 0; x(9) = 0;      % 17 stages, v=4
%     %x(1:9) = 0; x(11:14) = 0; x(16) = 0;   % 26 stages, v=5
%     %x(1:16) = 0; x(18:23) = 0; x(25) = 0;  % 37 stages, v=6
%     %x(1:25) = 0; x(27:34) = 0; x(36) = 0;  % 50 stages, v=7
%     %x(1:(v-2)^2) = 0; x((v-2)^2+2:end-5) = 0; x(end-4) = 0;
%     
%     %x(1:(v-2)^2) = 0; x((v-2)^2+2:end-4) = 0; % or this one
% 
%     %x(1:(v-2)^2) = 0; x((v-2)^2+2:end-3) = 0; % or this one
%     
%     
%     % plugging optimization variables in alpha matrix (for a full lower diagonal)
%     %alpha = diag(x(2:end-1),-2*v);
%     %alpha(2*v,1) = x(1);
%     %alpha(end,1) = 2/(s*((v-1)^2+1));
%     
%     % with only 2 variables
%     alpha = zeros(s+1,s+1);
%     %alpha(v^2-2*v+4,(v-2)^2) = x(1);
%     %alpha(end,v^2-2*v+2) = x(2);
%     %alpha(end,1) = 2/(s*((v-1)^2+1));
%     
%     % setting alpha(end,v^2-2*v+2) as free parameter
%     %r = -x(end); 
% %     alpha(v^2-2*v+4,(v-2)^2) = 1 - r*(v-1)/(s*(2*v-1)*(1-x(2)-2/(s*((v-1)^2+1))));
% %     x(1) = alpha(v^2-2*v+4,(v-2)^2);
% %     alpha(end,v^2-2*v+2) = x(2);
% %     alpha(end,1) = 2/(s*((v-1)^2+1));
%     
%     % setting alpha(v^2-2*v+4,(v-2)^2) as free parameter
% %     r = -x(end);
% %     alpha(v^2-2*v+4,(v-2)^2) = x(1);
% %     alpha(end-1,v^2-2*v+1) = 1 - 2/(s*((v-1)^2+1)) - r*(v-1)/(s*(2*v-1)*(1-x(1)));
% %     x(2) = alpha(end-1,v^2-2*v+1);
% %     alpha(end,1) = 2/(s*((v-1)^2+1));
% 
% 
% 
%     
%     alpha(v^2-2*v+4,(v-2)^2) = (v^2-1 - sqrt(v^3-3*v^2+v+1))/(4*v^2-6*v+2); %also (v^2-1 + sqrt(v^3-3*v^2+v+1))/(4*v^2-6*v+2)
%     alpha(end,1) = 2/(s*((v-1)^2+1));
%     alpha(end,v^2-2*v+2) = 1 - alpha(end,1) - (r*(v-1))/(s*(2*v-1)*(1-alpha(v^2-2*v+4,(v-2)^2)));
%     
%     alpha(end,1) = x(1);
%     alpha(end,v^2-2*v+2) = x(2);
%     
%     % Construction alpha matrix
%     big_diag = 1 - sum(alpha,2);
%     alpha = alpha + diag(big_diag(2:end),-1);
%     alpha = alpha(:,1:end-1); % fixing dimensions
%     
%     % Construct beta matrix
%     beta = alpha/r;
%     beta(end,1) = 0;
%     
% %     % Consistency for alpha matrix in the Shu-Osher form
% %     for i=2:s+1
% %        alpha(i,1)=1-sum(alpha(i,2:s));
% %     end
%     
% else


