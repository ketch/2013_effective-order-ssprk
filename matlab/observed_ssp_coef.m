function [sigma,r] = observed_ssp_coef(stages,effective_order,order, ...
    problem,initial_data,x_left,x_right,t_final,gridpts)
% Finds the observed SSP coefficient for a given problem
%
% Used by table_obs_SSP_coef.m

%==========================================================================

format long

% Set path
cd Source/

%==========================================================================

%% Editable parameters:

% For advection:
advection_speed = 1.0;

% For burgers 2nd order spatial discretization (burgers2):
    % available limiters:
        % 1. Koren
        % 2. Osher
        % 3. Superbee
lim_meth_burgers = 2;

% For buckley-leverett:
    % available limiters:
        % 1. Koren
        % 2. Osher
        % 3. Superbee
lim_meth_bl = 2;
alpha = 1./3;

% tolerence
tol = 10*eps;

%==========================================================================

%% Spatial discretization

% define dx
dx = (x_right - x_left)/gridpts;

% spacial mesh points
x = x_left:dx:x_right-dx;

% Define initial condition for y
if strcmp(initial_data,'cont')
    y_initial = 0.5 - 0.25*sin(pi.*x); 
elseif strcmp(initial_data,'discont')
    x_mid = (x_right + x_left)/2.;
    x_l = (x_left + x_mid)/2.;
    x_r = (x_right + x_mid)/2.;
    y_initial = 1*((x_l <= x) & (x <= x_r)) + 0 ; 
    if strcmp(problem,'buckley_leverett')
        x_mid = (x_right + x_left)/2.;
        y_initial = 0.5*((x >= x_mid)) + 0 ;
    end
end

% Set dt_FE and problem parameters
if strcmp(problem,'advection')
    dt_FE = dx/advection_speed;
    parameters = advection_speed;
elseif strcmp(problem,'burgers')
    dt_FE = dx/max(abs(y_initial));
    parameters = [];
elseif strcmp(problem,'burgers2')
    if lim_meth_burgers == 1
    	dt_FE = 0.006628140703518;
        if strcmp(initial_data,'discont')
            dt_FE = 0.004939698492462;
        end
    elseif lim_meth_burgers == 2
    	dt_FE = 0.008879396984925;
        if strcmp(initial_data,'discont')
            dt_FE = 0.006628140703518;
        end
    elseif lim_meth_burgers == 3
    	dt_FE = 0.006628140703518;
        if strcmp(initial_data,'discont')
            dt_FE = 0.004939698492462;
        end
    end
    parameters = lim_meth_burgers;
elseif strcmp(problem,'buckley_leverett')
    if (strcmp(initial_data,'discont') && gridpts == 100 && x_left == 0 ...
            && x_right == 1)
        if lim_meth_bl == 1
            dt_FE = 0.0024; % 0.002482412060302;
        elseif lim_meth_bl == 2
            dt_FE = 0.00304; % 0.003155778894472;
        elseif lim_meth_bl == 3
            dt_FE = 0.0024; % 0.002472361809045;
        end
        parameters = [alpha,lim_meth_bl];
    else
        disp('Correct parameters: discont,x_left=0,x_right=1,gridpts=100.')
        return;
    end
end

%==========================================================================

%% Temporal discretization

% Get butcher matrix
method = sprintf('../Methods/ESSPRK%d%d%d.mat',stages,effective_order, ...
    order);
load(method)

% Set auxiliary parameters
aux = [dx, parameters];

% define method's coefficient
sigma = r;

tvd_flag = 1;
while tvd_flag

    % Set time-step
    dt = sigma*dt_FE;
    
    % Compute initial TV-norm
    TV_initial = tvnorm(y_initial);

    % Starting method
    y = onestepRK(0,y_initial,dt,R,Rb,Rc,problem,aux);
    t = dt;
    
    % Check TVD condition
    tv_old = TV_initial;
    tv_new = tvnorm(y);
    if (tv_new > tv_old + tol)
        fprintf('ESSPRK%d%d%d:\n',stages,effective_order,order)
        fprintf('TV norm increase at %d \t tv_new - tv_old = %d\n', ... 
            t,tv_new - tv_old);
        sigma = sigma - 0.01;
        fprintf('Increase = %2.f%% \t\t\t\t coeff = %3.2f\n\n', ... 
            sigma/r*100-100,sigma)
        
        % Set path
        cd ..
        
        return;
    end
    tv_old = tv_new;

    while t <= t_final - dt
        % main method
        y = onestepRK(t,y,dt,M,b,c,problem,aux);
        %figure(2); clf; plot(x, y); drawnow()
        
        % Update time
        t = t + dt;
        
        % compute TVD norm
        tv_new = tvnorm(y);
        if (tv_new > tv_old + tol)
            fprintf('ESSPRK%d%d%d:\n',stages,effective_order,order)
            fprintf('TV norm increase at %d\t tv_new - tv_old = %d\n', ...
                t,tv_new - tv_old);
            sigma = sigma - 0.01;
            fprintf('Increase = %2.f%% \t\t\t\t coeff = %3.2f\n\n', ... 
            sigma/r*100-100,sigma)
            
            % Set path
            cd ..
        
            return;
        end
        tv_old = tv_new;

    end

    % stopping method
    %ysave = y;
    y = onestepRK(t,y,dt,T,Tb,Tc,problem,aux);
    
    % Final time
    t = t + dt;
    
    % Compute TV-norm for final solution
    tv_new = tvnorm(y);
    if (tv_new > tv_old + tol)
        fprintf('ESSPRK%d%d%d:\n',stages,effective_order,order)
        fprintf('TV norm increase at %d \t tv_new - tv_old = %d\n', ... 
            t,tv_new - tv_old);
        sigma = sigma - 0.01;
        fprintf('Increase = %2.f%% \t\t\t\t coeff = %3.2f\n\n', ... 
            sigma/r*100-100,sigma)
        
        % Set path
        cd ..
        
        return;
    end
    
    %figure(1); clf; plot(x, y); hold on; plot(x,ysave,'r.-');
    %t
    sigma = sigma + 0.01;
    %drawnow()
    %pause()
    
end

end
