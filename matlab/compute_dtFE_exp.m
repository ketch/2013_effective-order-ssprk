function [dtFE,r,dt] = compute_dtFE_exp(problem,initial_data,x_left, ...
    x_right,t_final,gridpts)
% Compute expected dtFE of TVD spatial discretizations with FE method 
% for various HCL problems
%
% Input variable meanings:
%    problem       --- 'advection','burgers','burgers2','buckley_leverett'
%    method        --- 'ESSPRKsqp', where s = stages, q = effective order, 
%                       p = order
%    initial_data  --- 'cont','discont', continuous and discontinuous data
%    x_left        --- left spatial point
%    x_right       --- right spatial point
%    t_final       --- final time
%    gridpts       --- number of spatial grid points
%
% Output variables:
%   dtFE           --- experimental SSP time-step for FE method
%   r              --- TV_norm ratio
%   dt             --- time-step range

%==========================================================================

format long

% Set path
cd Source/

%==========================================================================

%% Editable parameters:

% For advection:
advection_speed = 1.0;

% For burgers:
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
tol = 1e-14;

%==========================================================================

%% Spatial discretization

% define dx
dx = (x_right - x_left)/gridpts;

% spacial mesh points
x = x_left:dx:x_right-dx;

% Define initial condition for y
if strcmp(initial_data,'cont')
    y_initial = 1/2 - 1/4*sin(pi.*x);
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

% Set dt bounds and problem parameters
if strcmp(problem,'advection')
    parameters = advection_speed;
    dt_left = 0.005;
    dt_right = 0.015;
elseif (strcmp(problem,'burgers') || strcmp(problem,'burgers2'))
    parameters = lim_meth_burgers;
    dt_left = 0.001;
    dt_right = 0.015;
elseif strcmp(problem,'buckley_leverett')
    if (strcmp(initial_data,'discont') && gridpts == 100 && x_left == 0 ...
            && x_right == 1)
        parameters = [alpha,lim_meth_bl];
        dt_left = 0.002;
        dt_right = 0.004;
    else
        disp('Correct parameters: discont,x_left=0,x_right=1,gridpts=100.')
        return;
    end
end

%==========================================================================

%% Temporal discretization

% Forward Euler method
A = 0; c = 0; b = 1;

% Set time-step
dt = linspace(dt_left,dt_right,200);

% Set auxiliary parameters
aux = [dx, parameters];

%==========================================================================

%% Numerical solution

% initialize ratio
r = zeros(length(dt),1);
temp = 1;
for i = 2:length(dt)
    
    n = floor(t_final./dt(i));
    TV = zeros(n,1);
    
    y = y_initial;
    t = 0;
        
    % Compute initial TV-norm
    tv_old = tvnorm(y);

    for j=1:n
    
        % Main method
        y = onestepRK(t,y,dt(i),A,b,c,problem,aux);
    
        % Check TVD condition
        TV(j) = tvnorm(y)./tv_old;
        tv_old = tvnorm(y);
    
        % Update time
        t = t + dt(i);
    
    end
    
    r(i) = norm(TV,'inf');
     if (r(i) - 1.0 <= tol) && (i - temp == 1)
         dtFE = dt(i);
         temp = i;
     end

end

if dtFE == dt(end)
    disp('Need to increase dt right bound.')
end

%% Plot final solution:

fig = figure(1); clf;
plot(dt(2:end),r(2:end),'b-');
hc = get(fig,'children'); set(hc, 'fontsize', 14);
xlabel('$$\Delta t$$','FontSize',20,'Interpreter','latex');
ylabel('$$r$$','FontSize',20,'Interpreter','latex');

% Set path
cd ..

end