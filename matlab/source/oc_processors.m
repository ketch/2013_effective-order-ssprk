function coneq = oc_processors(x,Rs,Ts,M,b,c,effp,explicit)
% Order conditions for starting and stopping methods, R and T.
% The conditions are algebraic expressions on the Butcher coefficients
% and we follow Butcher's simplyfying result: beta_1 = 0
%
% Used by pre_post_processors.m and processors_nonlinear_con.m

%==========================================================================

% Extract coefficient for methods R,T and conditions on beta from x
if explicit == 1
    x1 = x(1:(Rs+1)*Rs/2);
    x2 = x(length(x1)+1:length(x1)+(Ts+1)*Ts/2);
else
    x1 = x(1:(Rs+3)*Rs/2);
    x2 = x(length(x1)+1:length(x1)+(Ts+3)*Ts/2);
end
x3 = x(length(x1)+length(x2)+1:end-3);

% Unpack methods R and T
[R,Rb] = unpack_RK(x1,Rs,explicit,0);
Rc = sum(R,2);
[T,Tb] = unpack_RK(x2,Ts,explicit,0);
Tc = sum(T,2);

% Pre-allogation of variables related to methods M, R, T and S
alpha = zeros(1,8);
rho = zeros(1,8);
tau = zeros(1,8);
beta = zeros(1,8);

%==========================================================================

%% Declarations for method M (alpha)

% order 1
alpha(1) = b'*ones(length(b),1);


% order 2
alpha(2) = b'*c;

% order 3
alpha(3) = b'*diag(c)*c;
alpha(4) = b'*M*c;

% order 4
alpha(5) = b'*diag(c)^2*c;
alpha(6) = b'*diag(M*c)*c;
alpha(7) = b'*M*diag(c)*c;
alpha(8) = b'*M^2*c;

% % order 5
% alpha(9) = b'*diag(c)^3*c;
% alpha(10) = b'*diag(c)^2*M*c;
% alpha(11) = b'*diag(c)*M*diag(c)*c;
% alpha(12) = b'*diag(c)*M^2*c;
% alpha(13) = b'*diag(M*c)*M*c;
% alpha(14) = b'*M*diag(c)^2*c;
% alpha(15) = b'*M*diag(c)*M*c;
% alpha(16) = b'*M^2*diag(c)*c;
% alpha(17) = b'*M^3*c;


%% Declarations for methods R and T (rho and tau respectively)

% order 1
rho(1) = Rb'*ones(Rs,1);                  tau(1) = Tb'*ones(Ts,1);


% order 2
rho(2) = Rb'*Rc;                          tau(2) = Tb'*Tc;

% order 3
rho(3) = Rb'*diag(Rc)*Rc;                 tau(3) = Tb'*diag(Tc)*Tc;
rho(4) = Rb'*R*Rc;                        tau(4) = Tb'*T*Tc;

% order 4
rho(5) = Rb'*diag(Rc)^2*Rc;               tau(5) = Tb'*diag(Tc)^2*Tc;
rho(6) = Rb'*diag(R*Rc)*Rc;               tau(6) = Tb'*diag(T*Tc)*Tc;
rho(7) = Rb'*R*diag(Rc)*Rc;               tau(7) = Tb'*T*diag(Tc)*Tc;
rho(8) = Rb'*R^2*Rc;                      tau(8) = Tb'*T^2*Tc;

% % order 5
% rho(9) = Rb'*diag(Rc)^3*Rc;             tau(9) = Tb'*diag(Tc)^3*Tc;
% rho(10) = Rb'*diag(Rc)^2*R*Rc;          tau(10) = Tb'*diag(Tc)^2*T*Tc;
% rho(11) = Rb'*diag(Rc)*R*diag(Rc)*Rc;   tau(11) = Tb'*diag(Tc)*T*diag(Tc)*Tc;
% rho(12) = Rb'*diag(Rc)*R^2*Rc;          tau(12) = Tb'*diag(Tc)*T^2*Tc;
% rho(13) = Rb'*diag(R*Rc)*R*Rc;          tau(13) = Tb'*diag(T*Tc)*T*Tc;
% rho(14) = Rb'*R*diag(Rc)^2*Rc;          tau(14) = Tb'*T*diag(Tc)^2*Tc;
% rho(15) = Rb'*R*diag(Rc)*R*Rc;          tau(15) = Tb'*T*diag(Tc)*T*Tc;
% rho(16) = Rb'*R^2*diag(Rc)*Rc;          tau(16) = Tb'*T^2*diag(Tc)*Tc;
% rho(17) = Rb'*R^3*Rc;                   tau(17) = Tb'*T^3*Tc;

%==========================================================================

%% Conditions on mehtod S (beta)

beta(1) = 0.0;
beta(2) = -1/6 + 1/2*alpha(3);

%% Effective order conditions on R and T

% Effective order 1
if effp >= 1
    % Conditions on R
    coneq(1) = rho(1) - alpha(1);                                          % rho(1) - (alpha(1) + beta(1));
    % Conditions on T
    coneq(2) = tau(1) - alpha(1);                                          % alpha(1) - (beta(1) + tau(1));
end

% Effective order 2
if effp >=2
    % Conditions on R                                      
    coneq(3) = rho(2) - (alpha(2) + beta(2));                              % rho(2) - (alpha(2) + beta(1)*alpha(1) + beta(2));
    % Conditions on T
    coneq(4) = tau(2) - (alpha(2) - beta(2));                              % alpha(2) - (beta(2) + tau(1)*beta(1) + tau(2));
end

% Effective order 3
if effp >= 3
    % Conditions on beta
    beta(3) = x3(1);
    beta(4) = x3(2);
    if effp > 3
        beta(3) = 1/12 - 1/2*alpha(3) + 1/3*alpha(5);
        beta(4) = -1/24 - 1/3*alpha(5) + alpha(6);
    end
    % Conditions on R
    coneq(5) = rho(3) - (alpha(3) + beta(3));                              % rho(3) - (alpha(3) + 2*beta(1)*alpha(2) + beta(1)^2*alpha(1) + beta(3));
    coneq(6) = rho(4) - (alpha(4) + alpha(1)*beta(2) + beta(4));           % rho(4) - (alpha(4) + beta(1)*alpha(2) + beta(2)*alpha(1) + beta(4));
    % Conditions on T
    coneq(7) = tau(3) - (alpha(3) - 2*alpha(1)*beta(2) - beta(3));         % alpha(3) - (beta(3) + 2*tau(1)*beta(2) + tau(1)^2*beta(1) + tau(3));
    coneq(8) = tau(4) - (alpha(4) - alpha(1)*beta(2) - beta(4));           % alpha(4) - (beta(4) + tau(1)*beta(2) + tau(2)*beta(1) + tau(4));
end

% Effective order 4
if effp >= 4
    % conditions on beta
    beta(5) = x3(1);
    beta(6) = x3(2);
    beta(7) = x3(3);
    beta(8) = x3(4);
    % Conditions on R
    coneq(9) = rho(5) - (alpha(5) + beta(5));                              % rho(5) - (alpha(5) +3*beta(1)*alpha(3) + 3*beta(1)^2*alpha(2) + beta(1)^3*alpha(1) + beta(5)); 
    coneq(10) = rho(6) - (alpha(6) + alpha(2)*beta(2) + beta(6));          % rho(6) - (alpha(6) + beta(1)*alpha(4) + beta(1)*alpha(3) + (beta(1)^2 + beta(2))*alpha(2) + beta(1)*beta(2)*alpha(1) + beta(6)); 
    coneq(11) = rho(7) - (alpha(7) + alpha(1)*beta(3) + beta(7));          % rho(7) - (alpha(7) + 2*beta(1)*alpha(4) + beta(1)^2*alpha(2) + beta(3)*alpha(1) + beta(7));
    coneq(12) = rho(8) - (alpha(8) + alpha(1)*beta(4) + ...
        alpha(2)*beta(2) + beta(8));                                       % rho(8) - (alpha(8) + beta(1)*alpha(4) + beta(2)*alpha(2) + beta(4)*alpha(1) + beta(8));
    % Conditions on T
    coneq(13) = tau(5) - (alpha(5) - 3*alpha(1)^2*beta(2) - ...
        3*alpha(1)*beta(3) - beta(5));                                     % alpha(5) - (beta(5) +3*tau(1)*beta(3) + 3*tau(1)^2*beta(2) + tau(1)^3*beta(1) + tau(5));
    coneq(14) = tau(6) - (alpha(6) - (alpha(1)^2 + alpha(2) - ... 
        beta(2))*beta(2) - alpha(1)*beta(3) - alpha(1)*beta(4) - beta(6)); % alpha(6) - (beta(6) + tau(1)*beta(4) + tau(1)*beta(3) + (tau(1)^2 + tau(2))*beta(2) + tau(1)*tau(2)*beta(1) + tau(6));
    coneq(15) = tau(7) - (alpha(7) - 2*alpha(1)*beta(4) - ...
        alpha(1)^2*beta(2) - beta(7));                                     % alpha(7) - (beta(7) + 2*tau(1)*beta(4) + tau(1)^2*beta(2) + tau(3)*beta(1) + tau(7));
    coneq(16) = tau(8) - (alpha(8) - alpha(1)*beta(4) - ...
        alpha(2)*beta(2) + beta(2)^2 - beta(8));                           % alpha(8) - (beta(8) + tau(1)*beta(4) + tau(2)*beta(2) + tau(4)*beta(1) + tau(8));
end

end
