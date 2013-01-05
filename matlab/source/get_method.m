function [R,Rb,Rc,M,b,c,T,Tb,Tc,r1,r,r2,effr1,ceff,effr2,...
    P,Pb,Pc,order_P] = get_method(stages,effective_order,order,class)
% Create ESSPRK scheme TMR, where R and T are SSP pre- and post-
% methods and M is the main ESSPRK(s,q,p) method.
% Also the "big" single method P = TMR is created using Butcher's 
% composition rule to verify the order is the desired one.

% TODO: Fix order conditions of implicit methods for orders more than 5
%==========================================================================

%% ESSPRK scheme

if nargin == 3
    % Creating main method of effective order q
    [M,b,c,r,ceff] = effective_ssp(stages,effective_order,order);
    % Creating starting and stopping methods R = SM and T = MS^(-1) for a 
    % given main method:
    [R,Rb,Rc,r1,effr1,T,Tb,Tc,r2,effr2] = pre_post_processors(M,b,c,r, ...
        effective_order,order);
else
    % Creating main method of effective order q
    [M,b,c,r,ceff] = effective_ssp(stages,effective_order,order,class);
    % Creating starting and stopping methods R = SM and T = MS^(-1) for a 
    % given main method:
    [R,Rb,Rc,r1,effr1,T,Tb,Tc,r2,effr2] = pre_post_processors(M,b,c,r, ...
        effective_order,order,class);
end

%==========================================================================

%% P = TMR method

% # of stages of consisting methods
s1 = length(Rb);
s2 = length(b);
s3 = length(Tb);

% Contructing P = R*M*T to verify the classical order of method P
P = [R zeros(s1,s2) zeros(s1,s3) ; ones(s2,1)*Rb' M zeros(s2,s3) ; ... 
    ones(s3,1)*[Rb ; b]' T];
Pb = [Rb; b ; Tb];
Pc = sum(P,2);

%% Order conditions

% order 1
m(1) = Pb'*ones(length(Pb),1);
% order 2
m(2) = Pb'*Pc;
% order 3
m(3) = Pb'*diag(Pc)*Pc;
m(4) = Pb'*P*Pc;
% order 4
m(5) = Pb'*diag(Pc)^2*Pc;
m(6) = Pb'*diag(P*Pc)*Pc;
m(7) = Pb'*P*diag(Pc)*Pc;
m(8) = Pb'*P^2*Pc;
% order 5
m(9) = Pb'*diag(Pc)^3*Pc;
m(10) = Pb'*diag(Pc)^2*P*Pc;
m(11) = Pb'*diag(Pc)*P*diag(Pc)*Pc;
m(12) = Pb'*diag(Pc)*P^2*Pc;
m(13) = Pb'*diag(P*Pc)*P*Pc;
m(14) = Pb'*P*diag(Pc)^2*Pc;
m(15) = Pb'*P*diag(Pc)*P*Pc;
m(16) = Pb'*P^2*diag(Pc)*Pc;
m(17) = Pb'*P^3*Pc;

%% Checking order of method P

% Order 2
if effective_order == 2
    % B-series of E^3
    e = [3. 9./2];
    if norm(m(1:2) - e,'inf') < 1e-14
        order_P = 2;
    else
        order_P = 'Not order 2';
    end
end

% Order 3
if effective_order == 3
    % B-series of E^3
    e = [3. 9./2 9. 9./2];
    if norm(m(1:4) - e,'inf') < 1e-14
        order_P = 3;
    else
        order_P = 'Not order 3';
    end
end

% Order 4
if effective_order == 4
    % B-series of E^3
    e = [3. 9./2 9. 9./2 81./4 81./8 81./12 81./24];
    if norm(m(1:8) - e,'inf') < 1e-14
        order_P = 4;
    else
        order_P = 'Not order 4';
    end
end

end