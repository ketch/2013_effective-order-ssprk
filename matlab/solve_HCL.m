function solve_HCL(problem,method,initial_data,x_left,x_right,t_final, ...
    gridpts,sav,movie)
% Solve advection equation with first order upwind spatial disretization
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
%    sav           --- 1 to save to file, 0 not
%    movie         --- 1 to creat movie, 0 not

%==========================================================================

format long

% Set path
cd Source/

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


%% Spatial discretization

% define dx
dx = (x_right - x_left)/gridpts;

% spacial mesh points
x = x_left:dx:x_right-dx;

% Define initial condition for y
if strcmp(initial_data,'cont')
    y_initial = 0.5 - 0.25*sin(pi.*x);
    setaxis = [x_left x_right 0.23 0.78];
    %setaxis = [.38 .82 .645 .745];
elseif strcmp(initial_data,'discont')
    x_mid = (x_right + x_left)/2.;
    x_l = (x_left + x_mid)/2.;
    x_r = (x_right + x_mid)/2.;
    y_initial = 1*((x_l <= x) & (x <= x_r)); 
    setaxis = [x_left x_right -0.1 1.1];
    %setaxis = [1.25 1.8 0.997 1.0004];
    if strcmp(problem,'buckley_leverett')
        x_mid = (x_right + x_left)/2.;
        y_initial = 0.5*((x >= x_mid));
        setaxis = [x_left x_right -0.1 0.6];
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
    
%% Temporal discretization

% Get butcher matrix
method = sprintf('../Methods/%s',method);
load(method);

% Set SSP coefficient
coef = r;

% Set time-step
dt = coef*dt_FE;

%% Numerical solution

if movie == 1
    mov = figure(1);
    winsize = get(mov,'Position');
    winsize(1:2) = [0 0];
    numframes = 150;
    A = moviein(numframes,mov,winsize);
    set(mov,'NextPlot','replacechildren')
    i=1;
end

% Set auxiliary parameters
aux = [dx, parameters];

% Compute initial TV-norm
tv_old = tvnorm(y_initial);

% Starting method
y = onestepRK(0,y_initial,dt,R,Rb,Rc,problem,aux);

% Update time
t = dt;

% Check TVD condition
tv_new = tvnorm(y);
if (tv_new > tv_old + tol)
    fprintf('TV norm increase at t = dt = %d\n', t);
    fprintf('tv_new - tv_old = %d\n', tv_new - tv_old);
    fprintf('tv_old = %d\n', tv_old)
    fprintf('tv_new = %d\n\n', tv_new)
    pause;
end
tv_old = tv_new;

% Main computation
while t < t_final - dt
%for n=2:floor(t_final/dt)-1
    
    % Main method
    y = onestepRK(t,y,dt,M,b,c,problem,aux);
    
    % Update time
    t = t + dt;
    
    % newx = 0*(x < 0.5) + ((x-0.5)/t).*((0.5 <= x) & (x <= 0.5 + t)) ...
    % + 1*((0.5 + t < x) & (x <= 1.5 + t/2)) + 0*(x > 1.5 + t/2);
    % newx = x + advection_speed*t; % actual solution by method of ... 
    % characteristics for advection
    % newx = x + y_initial*t; % actual solution by method of ... 
    % characteristics for burgers
    
    % Compute TV-norm
    tv_new = tvnorm(y);

    % Plot solution
    h = figure(1); clf;
    plot(x, y,'r.-');
    hc = get(h,'children'); set(hc, 'FontSize', 14);
    hold on;
    title(['t=' num2str(t), '$$\quad$$ TV=' num2str(tv_new), ... 
        '$$\quad$$ TV difference=' num2str(tv_new - tv_old)], ...
        'FontSize',16,'Interpreter','latex');
    axis(setaxis);
    xlabel('$$x$$','FontSize',20,'Interpreter','latex');
    ylabel('$$u$$','FontSize',20,'Interpreter','latex');axis(setaxis);
    %h = plot(x,newx,'b--'); % for disc solution newx above
    %h = plot(newx,y_initial,'b--');
    %if norm(newx-sort(newx)) ~= 0, set(h,'color','g'), end
    
    % Check TVD condition
    if (tv_new > tv_old + tol)
        fprintf('TV norm increase at t = %d\n', t);
        fprintf('tv_new - tv_old = %d\n', tv_new - tv_old);
        fprintf('tv_old = %d\n', tv_old);
        fprintf('tv_new = %d\n\n', tv_new);
        pause;
    end
    tv_old = tv_new;
    
    if movie ==1
        if i <= numframes;
            A(:,i) = getframe(mov,winsize); 
            i= i+1;
        end
    end
    
end

% Stopping method
y = onestepRK(t,y,dt,T,Tb,Tc,problem,aux);

% Final time
t = t + dt;

% Compute TV-norm for final solution
tv_new = tvnorm(y);

% Check TVD condition
if (tv_new > tv_old + tol)
    fprintf('TV norm increase at t_final = %d\n', t);
    fprintf('tv_new - tv_old = %d\n', tv_new - tv_old);
    fprintf('tv_old = %d\n', tv_old)
    fprintf('tv_new = %d\n\n', tv_new)
    pause;
end

% Plot initial function of y 
hold on;
plot(x, y_initial);
legend('numerical solution','initial solution')


%% Save movie

if movie == 1
    file = sprintf('../Movies/%s_%s.avi',problem, initial_data);
    movie2avi(A, file, 'compression', 'none', 'fps', 6);
end


%% Plot final solution:
%setaxis = [.95 1.65 0.997 1.0005];
fig = figure(2); clf;
plot(x, y,'r*-');
hc = get(fig,'children'); set(hc, 'fontsize', 14);
title(['$$t =\; $$' num2str(t),'$$\qquad TV =\; $$' ... 
    num2str(tvnorm(y))],'FontSize',16,'Interpreter','latex');
axis(setaxis);
xlabel('$$x$$','FontSize',20,'Interpreter','latex');
ylabel('$$u$$','FontSize',20,'Interpreter','latex');

%% Save to file

if sav == 1
    file = sprintf('../Figures/%s_%s.eps',problem, initial_data);
    print('-depsc', file)
end

% Set path
cd ..