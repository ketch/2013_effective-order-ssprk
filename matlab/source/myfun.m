function F = myfun(t,y,problem,aux)
% Several RHS for ODE problems
%
% Used by onestepRK.m

%==========================================================================

%% Define RHS 

if strcmp(problem,'linear')
    lambda = aux; F = lambda.*y;
elseif strcmp(problem,'nonlinear1')
    lambda = aux(1); r = aux(2); F = lambda.*y.^r;
elseif strcmp(problem,'nonlinear2')
    F = -y.*(1-y);
elseif strcmp(problem,'nonlinear3')
    d = aux; F = (cos(d.*y)).^2;
elseif strcmp(problem,'nonlinear4')
    F = 4.*y.*(sin(t)).^3.*cos(t);
elseif strcmp(problem,'nonlinear5')
    F = 4.*t.*sqrt(y);
elseif strcmp(problem,'nonlinear6')
    F = y./t.*log(y);
elseif strcmp(problem,'nonlinear7')
    F = nonlinear7(t,y);
elseif strcmp(problem,'ode1')
    F = ode1(t,y);
    F = F';
elseif strcmp(problem,'ode2')
    F = ode2(t,y);
    F = F';
elseif strcmp(problem,'ode3')
    F = ode3(t,y);
    F = F';
elseif strcmp(problem,'vdp')
    F = vdp(t,y);
    F = F';
elseif strcmp(problem,'krogh')
    F = krogh(t,y);
    F = F';
elseif strcmp(problem,'stanescu1998')
    F = stanescu1998(t,y);
    F = F';
elseif strcmp(problem,'advection')
    % parameters
    a = aux(2);
    
    % grid-spacing and initialization
    dx = aux(1);
    
    % flux
    f = @(u,t) a.*u;
    
    % periodic boundary conditions
    % U_{0}=y(m),U_{1}=y(1),...,U_{m-1}=y(m-1),U_{m}=y(m)
    uj = y(1:end); % u_{j}, j = 1,...,m
    ujm1 = y([end,1:end-1]); % u_{j-1}, j = 1,...,m

    % numerical flux
    F = (-1./dx)*(f(uj) - f(ujm1)); % F = -1/(dx).*(uj - ujm1);
elseif strcmp(problem,'burgers')
    % grid-spacing and initialization
    dx = aux(1);

    % flux
    f = @(u,t) 0.5*u.^2;
    
    % periodic boundary conditions
    % U_{0}=y(m),U_{1}=y(1),...,U_{m-1}=y(m-1),U_{m}=y(m)
    uj = y(1:end); % u_{j}, j = 1,...,m
    ujm1 = y([end,1:end-1]); % u_{j-1}, j = 1,...,m
    
    % numerical flux
    F = (-1./dx)*(f(uj) - f(ujm1)); % F(2:end) = -1/(2*dx)*(uj.^2-ujm1.^2);
elseif strcmp(problem,'burgers2')
    % grid-spacing and initialization
    dx = aux(1);
    limiter_method = aux(2);
    
    % flux
    f = @(u,t) 0.5*u.^2;
    
    % periodic boundary conditions
    % U_{0}=y(m),U_{1}=y(1),...,U_{m-1}=y(m-1),U_{m}=y(m)
    ujp1 = y([2:end,1]);       % u_{j+1}, j = 1,...,m
    uj = y(1:end);             % u_{j}, j = 1,...,m
    ujm1 = y([end,1:end-1]); % u_{j-1}, j = 1,...,m
    
    % limiter argument
    % theta = (uj - ujm1)/(ujp1 - uj), j = 1,...,m
    theta = zeros(size(uj));
    for i=1:length(uj);
       if(ujp1(i) - uj(i) == 0)
 	theta(i) = 1e6;
       else
 	theta(i) = (uj(i) - ujm1(i))./(ujp1(i) - uj(i));
       end
    end
    
    if limiter_method == 1
        % Koren limiter
        limiter = max(0, min(2, min(2./3. + 1./3.*theta, 2*theta)));
    elseif limiter_method == 2
        % Osher limiter
        beta = 1.5;
        limiter = max(0, min(theta,beta));
    elseif limiter_method == 3
        % Superbee limiter
        temp1 = max(0, min(theta,2));
        temp2 = max(0, min(2*theta,1));
        limiter = max(temp1, temp2);
    end
    
    % Numerical flux arguments
    ujph = uj + 0.5*limiter.*(ujp1 - uj);  % u_{j+1/2}, j = 1,...,m
    ujmh = ujph([end,1:end-1]);            % u_{j-1/2}, j = 1,...,m
    
    % numerical flux
    F = (-1./dx)*(f(ujph) - f(ujmh));
elseif strcmp(problem,'buckley_leverett')
    % parameters
    alpha = aux(2);
    limiter_method = aux(3);
    
    % grid-spacing and initialization
    dx = aux(1);
    
    % flux
    f = @(u,t) u.^2./(u.^2 + alpha*(1-u).^2);
    
    % periodic boundary conditions
    % U_{0}=y(m),U_{1}=y(1),...,U_{m-1}=y(m-1),U_{m}=y(m)
    ujp1 = y([2:end,1]);       % u_{j+1}, j = 1,...,m
    uj = y(1:end);             % u_{j}, j = 1,...,m
    ujm1 = y([end,1:end-1]); % u_{j-1}, j = 1,...,m
    
    % limiter argument
    % theta = (uj - ujm1)/(ujp1 - uj), j = 1,...,m
    theta = zeros(size(uj));
    for i=1:length(uj);
       if(ujp1(i) - uj(i) == 0)
 	theta(i) = 1e6;
       else
 	theta(i) = (uj(i) - ujm1(i))./(ujp1(i) - uj(i));
       end
    end
    
    if limiter_method == 1
        % Koren limiter
        limiter = max(0, min(2, min(2./3. + 1./3.*theta, 2*theta)));
    elseif limiter_method == 2
        % Osher limiter
        beta = 1.5;
        limiter = max(0, min(theta,beta));
    elseif limiter_method == 3
        % Superbee limiter
        temp1 = max(0, min(theta,2));
        temp2 = max(0, min(2*theta,1));
        limiter = max(temp1, temp2);
    end
    
    % Numerical flux arguments
    ujph = uj + 0.5*limiter.*(ujp1 - uj);  % u_{j+1/2}, j = 1,...,m
    ujmh = ujph([end,1:end-1]);            % u_{j-1/2}, j = 1,...,m
    
    % numerical flux
    F = (-1./dx)*(f(ujph) - f(ujmh));
end

end
