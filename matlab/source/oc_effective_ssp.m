function coneq = oc_effective_ssp(x,s,p,effp,explicit,shuoshercoefficients)
% Effective order conditions for Runge--Kutta methods.
% The conditions are algebraic expressions on the Butcher coefficients
% and we follow Butcher's simplyfying result: beta_1 = 0
%
% Used by effective_ssp.m and ssp_nonlinear_con.m

%==========================================================================

if shuoshercoefficients == 0
    % Extract the Butcher coefficients from x
    [A,b] = unpack_RK(x,s,explicit,shuoshercoefficients);
    c = sum(A,2);
else
    % Extract the Shu-0sher coefficients from x
    [alpha,beta] = unpack_RK(x,s,explicit,shuoshercoefficients);
    
    % Construct Butcher tableau
    [A,b] = shuosher2butcher(alpha,beta);
    c = sum(A,2);
end

% Pre-allogation
s = length(b); % # of stages
unity = ones(s,1);
m = zeros(17,1);

%==========================================================================

%% Declarations for method M

% order 1
m(1) = b'*unity;

% order 2
m(2) = b'*c;

% order 3
m(3) = b'*diag(c)*c;
m(4) = b'*A*c;

% order 4
m(5) = b'*diag(c)^2*c;
m(6) = b'*diag(A*c)*c;
m(7) = b'*A*diag(c)*c;
m(8) = b'*A^2*c;

% order 5
m(9) = b'*diag(c)^3*c;
m(10) = b'*diag(c)^2*A*c;
m(11) = b'*diag(c)*A*diag(c)*c;
m(12) = b'*diag(c)*A^2*c;
m(13) = b'*diag(A*c)*A*c;
m(14) = b'*A*diag(c)^2*c;
m(15) = b'*A*diag(c)*A*c;
m(16) = b'*A^2*diag(c)*c;
m(17) = b'*A^3*c;

%==========================================================================

%% Effective order 1
if effp >= 1
    coneq(1) = m(1) - 1;
end

%% Effective order 2
if effp >= 2 && p >= 2
    coneq(2) = m(2) - 1/2;
end

%% Effective order 3
if effp >= 3 && p >= 2
    
    % classical order 2 
    coneq(3) = m(4) - 1/6;
    
    if p >= 3
        % classical order 3 
        coneq(4) = m(3) - 1/3;
    end
    
end

%% Effective order 4
if effp >= 4 && p >= 2
    
    % classical order 2 to 4 
    coneq(5) = m(8) - 1/24;
    coneq(6) = - m(3) + m(5) - 2*m(6) + m(7) + 1/4;
    
    if p >= 4
        % classical order 4 
        coneq(6) = m(5) - 1/4;
        coneq(7) = m(6) - 1/8;
        coneq(8) = m(7) - 1/12;
    end
    
end

%% Effective order 5
if effp >= 5 && p >= 2
    
    beta_2 = 1/2*m(3) - 1/6;
    
    % classical order 2 to 5 
    coneq(9) = m(17) - 1/120;
    coneq(10) = -(beta_2)^2 + 1/4*m(9) - m(10) + m(13);
    coneq(11) = 6*beta_2^2 + 3/2*m(3) - m(5) - 1/2*m(9) + 3*m(10) - ...
        3*m(11) + m(14) - 3/10;
    coneq(12) = 2*beta_2^2 - m(11) - m(12) + m(15) + 1/2*m(3) - m(6) - ...
        1/2*m(9) + 2*m(10) - 1/15;
    coneq(13) = - 4*beta_2^2 - m(3) + m(5) - 2*m(6) + m(11) - 2*m(12) + ...
        m(16) + 19/60;
    
    if p >= 5
        % classical order 5 
        coneq(10) = m(9) - 1/5;
        coneq(11) = m(10) - 1/10;
        coneq(12) = m(11) - 1/15;
        coneq(13) = m(12) - 1/30;
        coneq(14) = m(13) - 1/20;
        coneq(15) = m(14) - 1/20;
        coneq(16) = m(15) - 1/40;
        coneq(17) = m(16) - 1/60;
    end
    
end

end
