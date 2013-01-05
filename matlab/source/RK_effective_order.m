function effp = RK_effective_order(A,b,c)
% Determines effective order of a RK method given the classical order,
% up to fifth order

% For an s-stage method, input A should be a s x s matrix; 
% b and c should are column vectors of length s

%==========================================================================

tol = 1.e-15;
s = length(b); % number of stages
unity = ones(s,1);
effp = 0;

%==========================================================================

%% Declarations for method A

% order 1
a(1) = b'*unity;

% order 2
a(2) = b'*c;

% order 3
a(3) = b'*diag(c)*c;
a(4) = b'*A*c;

% order 4
a(5) = b'*diag(c)^2*c;
a(6) = b'*diag(A*c)*c;
a(7) = b'*A*diag(c)*c;
a(8) = b'*A^2*c;

% order 5
a(9) = b'*diag(c)^3*c;
a(10) = b'*diag(c)^2*A*c;
a(11) = b'*diag(c)*A*diag(c)*c;
a(12) = b'*diag(c)*A^2*c;
a(13) = b'*diag(A*c)*A*c;
a(14) = b'*A*diag(c)^2*c;
a(15) = b'*A*diag(c)*A*c;
a(16) = b'*A^2*diag(c)*c;
a(17) = b'*A^3*c;

%==========================================================================

%% Effective order conditions

% Order 1
t = a(1) - 1;

if max(abs(t)) < tol
    effp = 1;
end

% Order 2
t = a(2) - 1/2;

if max(abs(t)) < tol && effp == 1
    effp = 2;
end

% Order 3   
t = a(4)- 1/6;

if max(abs(t)) < tol && effp == 2
    effp = 3;
end

% Order 4
t(1) = - a(3) + a(5) - 2*a(6) + a(7) + 1/4;
t(2) = a(8) - 1/24;

if max(abs(t)) < tol && effp == 3
    effp = 4;
end

% Order 5
beta_2 = a(3)/2 - 1/6;
t(1) = -(beta_2)^2 + 1/4*a(9) - a(10) + a(13);
t(2) = 6*beta_2^2 + 3/2*a(3) - a(5) - 1/2*a(9) + 3*a(10) - 3*a(11) + ...
    a(14) - 3/10;
t(3) = 2*beta_2^2 - a(11) - a(12) + a(15) + 1/2*a(3) - a(6) - ...
    1/2*a(9) + 2*a(10) - 1/15;
t(4) = - 4*beta_2^2 - a(3) + a(5) - 2*a(6) + a(11) - 2*a(12) + ...
    a(16) + 19/60;
t(5) = a(17) - 1/120;

if max(abs(t)) < tol && effp == 4
    effp = 5;
    disp('This method has at least effective order 5');
end

end

