function [err,h] = test_method(stages,effective_order,order,problem,plot,m)
% Solve initial value problem to test the order of an ESSPRK scheme.
% Problem options span from scalar linear to nonlinear scalar and system
% problems as described below.
%
% Input variable meanings:
%   stages          --- method's stages
%   effective_order --- method's effective order
%   order           --- method's classical order 
%   problem         --- 'advection','burgers','burgers2','buckley_leverett'
%   initial_data    --- 'cont','discont', continuous and discontinuous data
%   plot            --- 1 to creat plot, 0 not
%   m               --- number of spatial points
%
% Output variables:
%   err             --- error
%   h               --- time-step
%
% Detailed description of available problems:
%
% 1. Linear scalar
% dy/dt=lambda*y for t in (0,1), u(0)=1                      --> linear
%
% 2. Nonlinear scalar
% dy/dt=lambda*y^k for t in [0,1], u(0)=((1-k)*c)^(1/(1-k))  --> nonlinear1
% dy/dt=-y*(1-y) for t in [0,1], u(0)=1/(1-A);               --> nonlinear2
% dy/dt=(cos(d*y))^2 for t in [0,1], u(0)=atan(d*c)/d;       --> nonlinear3
% dy/dt=4*y*(sin(t))^3*cos(t) for t in [0,1], u(0)=1         --> nonlinear4
% dy/dt=4*t*sqrt(y) for t in [0,5], u(0)=1                   --> nonlinear5
% dy/dt=y/t*log(y) for t in [0.5,1], u(0)=exp(1)             --> nonlinear6
% dy/dt=(y-t)/(y+t) for t in [0,20], u(0)=5                  --> nonlinear7
%
% 3. Nonlinear systems
% dydt=[y(2);a*(1-y(1)^2)*y(2)-y(1)] for t in [0,50], y(0)=2, y'(0)=1 
%                                                            --> vdp
% dydt=[2*(y(1)-y(1)*y(2));-(y(2)-y(1)*y(2))]for t in [0,20],y(0)=1,y'(0)=3 
%                                                            --> ode1
% dydt=[y(1);2*y(2)] for t in [0,5], y(0)=1, y'(0)=1         --> ode2
% dydt=[-y(1);y(1)-y(2)^2;y(2)^2] for t in [0,20], y(0)=1,y'(0)=0,y''(0)=0
%                                                            --> ode3
% dydt=[y(2)*y(3);-y(1)*y(3);-0.51*y(1)*y(2)] for t in [0,12], 
% y(0)=0,y'(0)=1,y''(0)=1                                    --> krogh
% dydt=[1/y(1)-y(2)*exp(t^2)/t^2-t;1/y(2)-exp(t^2)-2*t*exp(-t^2)]
% for t in [1,1.4], y(0)=1, y'(0)=exp(-1)                  --> stanescu1998
%
% Used by convergence_test.m

%==========================================================================

format long

% Set path
cd Source/

% number of spacial points
if nargin == 5
    m = 2.^(1:6)*50;
end

%==========================================================================

%% Set method
method = sprintf('../Methods/ESSPRK%d%d%d.mat',stages,effective_order, ...
    order);
load(method)

%==========================================================================

%% Choose problem parameters

if (strcmp(problem,'advection') || strcmp(problem,'burgers') || ... 
        strcmp(problem,'buckley_leverett'))
    disp('No convergence study for PDEs.')
    return;
end
options = odeset('RelTol',1e-13,'AbsTol',1e-13*ones(1, ... 
    length(@(dimensions) dimensions)));

% Define exact (or reference) solution and initial conditions
if strcmp(problem,'linear')
    dimensions = 1;
    t_start = 0.;
    t_final = 1.;
    lambda = -1.5; exact = @(t) exp(lambda.*t);
    y_initial = 1;
    aux = lambda;
elseif strcmp(problem,'nonlinear1')
    dimensions = 1;
    t_start = 0.;
    t_final = 1.;
    lambda = -1.5; d= -0.1; r = 2; 
    exact = @(t) ((1-r)*(lambda.*t+d)).^(1/(1-r));
    y_initial = ((1-r)*d)^(1/(1-r));
    aux = [lambda,r];
elseif strcmp(problem,'nonlinear2')
    dimensions = 1;
    t_start = 0.;
    t_final = 1.;
    A = -100; exact = @(t) 1./(1-A.*exp(t));
    y_initial = 1/(1-A);
    aux = 'none';
elseif strcmp(problem,'nonlinear3')
    dimensions = 1;
    t_start = 0.;
    t_final = 1.;
    d1 = -1.5; d2 = 2.; exact = @(t) atan(d2.*(t+d1))./d2;
    y_initial = atan(d1*d2)/d2;
    aux = d2;
elseif strcmp(problem,'nonlinear4')
    dimensions = 1;
    t_start = 0.;
    t_final = 5.;
    u0 = 1.; exact = @(t) u0.*exp((sin(t)).^4);
    y_initial = u0;
    aux = 'none';
elseif strcmp(problem,'nonlinear5')
    dimensions = 1;
    t_start = 0.;
    t_final = 5.;
    exact = @(t) (1. + t.^2).^2;
    y_initial = 1.;
    aux = 'none';
elseif strcmp(problem,'nonlinear6')
    dimensions = 1;
    t_start = 0.5;
    t_final = 1.;
    exact = @(t) exp(2.*t);
    y_initial = exp(1.);
    aux = 'none';
elseif strcmp(problem,'nonlinear7')
    dimensions = 1;
    t_start = 0.0;
    t_final = 20.;
    y_initial = 5.;
    [t,y] = ode45('nonlinear7',[0 t_final],y_initial,options);
    exact = y;
    aux = 'none';
elseif strcmp(problem,'ode1')
    % Launch Matlab's ode45 routine
    % to solve ode system equation on [0,20] 
    % with initial condition y(0)=1, y'(0)=3
    dimensions = 2;
    t_start = 0.;
    t_final = 20;
    y_initial = [1. 3.]';
    [t,y] = ode45('ode1',[0 t_final],y_initial,options);
	exact = y(:,1);
    aux = 'none';
elseif strcmp(problem,'ode2')
    dimensions = 2;
    t_start = 0.;
    t_final = 5;
    y_initial = [1. 1.]';
    u = {@(t) exp(t), @(t) exp(2.*t)};
    exact = u{1};
    aux = 'none';
elseif strcmp(problem,'ode3')
    % Launch Matlab's ode45 routine
    % to solve ode system equation on [0,20] 
    % with initial condition y(0)=1, y'(0)=0, y''(0)=0
    dimensions = 3;
    t_start = 0.;
    t_final = 20;
    y_initial = [1.,0.,0.]';
    [t,y] = ode45('ode3',[0 t_final],y_initial,options);
	exact = y(:,1);
    aux = 'none';
elseif strcmp(problem,'vdp')
    % Launch Matlab's ode45 predictor-corrector routine
    % to solve van der Pol's equation on [0,50] 
    % with initial condition y(0)=2, y'(0)=1
    dimensions = 2;
    t_start = 0.;
    t_final = 50.;
    y_initial = [2. 1.]';
    [t,y] = ode45('vdp',[0 t_final],y_initial,options);
	exact = y(:,1);
    aux = 'none';
elseif strcmp(problem,'krogh')
    % Launch Matlab's ode45 routine
    % to solve Euler's ridig body equations on [0,12] 
    % with initial condition y(0)=0, y'(0)=1, y''(0)=1
    dimensions = 3;
    t_start = 0.;
    t_final = 12.;
    y_initial = [0. 1. 1.]';
    [t,y] = ode45('krogh',[0 t_final],y_initial,options);
	exact = y(:,1);
    aux = 'none';
elseif strcmp(problem,'stanescu1998')
    dimensions = 2;
    t_start = 1.;
    t_final = 1.4;
    y_initial = [1. exp(-1.)]';
    u = {@(t) 1./t, @(t) exp(-t.^2)};
    exact = u{1};
    aux = 'none';
end

% Define initial conditions for t
left_t = t_start;
right_t = t_final;

%==========================================================================

%% Compute numerical solution

% Set path
cd ../Source/

% Preallocations
h = zeros(length(m),1);
err = zeros(length(m),1);

for n = 1:length(m)

% Preallocation of variables
y = zeros(m(n),dimensions);

% Define mesh points
t = linspace(left_t,right_t,m(n));

% Initial conditions
y(1,:) = y_initial;

% Define step
h(n) = t(2) - t(1);

%% Numerical solution

% Starting method
y(2,:) = onestepRK(t(1),y(1,:),h(n),R,Rb,Rc,problem,aux);

% Main method
for k = 3:m(n)-1
    y(k,:) = onestepRK(t(k-1),y(k-1,:),h(n),M,b,c,problem,aux);
end

% Finishing method
y(end,:) = onestepRK(t(end-1),y(end-1,:),h(n),T,Tb,Tc,problem,aux);

% error estimation
if ~strcmp(problem,'nonlinear7') && (dimensions == 1 ... 
        || strcmp(problem,'ode2') || strcmp(problem,'stanescu1998')) 
    err(n) = abs(exact(t(end)) - y(end,1));
else
    err(n) = abs(exact(end) - y(end,1));
end

end

%==========================================================================

%% Plot of error
if plot == 1
    fig = figure(4); clf;
    loglog(h, h.^1, 'g--', h, h.^2, 'b--', h, h.^3, 'm--', ... 
        h, h.^4, 'r--', h, h.^5, 'c--', h, err, 'k-*');
    method = sprintf('Error when M is SSPRK(%d,%d,%d)', ... 
        stages,effective_order,order);
    legend('Order 1', 'Order 2', 'Order 3', 'Order 4', 'Order 5', ... 
        method, 'Location', 'SouthEast');
    xlabel('\Delta t');
    ylabel('Error')
    title('RM^{n-2}T')
    hold off
    pause;
    close(fig)
end

% Set path
cd ..

end