function convergence_test(problem,effective_order,m)
% Convergence test on linear and nonlinear ODEs for ESSPRK methods.
% Results are shown in a figure.
%
% Input variable meanings:
%   problem         --- linear and nonlinear problems - check test_method.m
%                       for more details
%   effective_order --- method's effective order
%   m               --- number of spatial points

%==========================================================================

if nargin == 2
    m = 2.^(1:6)*50;
end

if strcmp(problem,'vdp')
    m = 2.^(2:7)*100;
end

if effective_order == 3
    stages1 = 3; stages2 = 5; stages3 = 7; stages4 = 9;
    err1 = test_method(stages1,3,2,problem,0,m);
    err2 = test_method(stages2,3,2,problem,0,m);
    err3 = test_method(stages3,3,2,problem,0,m);
    [err4,h] = test_method(stages4,3,2,problem,0,m);
    order = '3rd order';
    stages1 = sprintf('$s=%d$',stages1);
    stages2 = sprintf('$s=%d$',stages2);
    stages3 = sprintf('$s=%d$',stages3);
    stages4 = sprintf('$s=%d$',stages4);
elseif effective_order == 4
    stages1 = 4; stages2 = 6; stages3 = 8; stages4 = 10;
    err1 = test_method(stages1,4,2,problem,0,m);
    err2 = test_method(stages2,4,2,problem,0,m);
    err3 = test_method(stages3,4,2,problem,0,m);
    [err4,h] = test_method(stages4,4,2,problem,0,m);
    order = '4th order';
    stages1 = sprintf('$s=%d$',stages1);
    stages2 = sprintf('$s=%d$',stages2);
    stages3 = sprintf('$s=%d$',stages3);
    stages4 = sprintf('$s=%d$',stages4);
end

save('Data/error.mat','err1','err2','err3','err4','h')

% Plot of errors
fig = figure(1); clf;
p = loglog(h, (2*h).^effective_order, 'r--', h, err1, 'b*-', h, err2, ...
    'bo-', h, err3,'bx-', h, err4,'bs-');
set(p, 'MarkerSize', 8);
hc = get(fig,'children'); set(hc, 'FontSize', 14);
leg = legend('String',{order,stages1,stages2,stages3,stages4});
set(leg,'Location','SouthEast','FontSize',20,'Interpreter','latex');
h = get(leg,'Position'); h(3) = 0.02+h(3); h(1) = -0.02+h(1);
set(leg,'Position',h);
xlabel('$\Delta t$','FontSize',20,'Interpreter','latex');
ylabel('$|error|$','FontSize',20,'Interpreter','latex');
hold on

file = sprintf('Figures/convergence_%d_order.eps',effective_order);
print('-depsc', file)

end

