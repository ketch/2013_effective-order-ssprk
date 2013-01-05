function p = RK_order(A,b,c)
% Determines order of a RK method, up to sixth order

% For an s-stage method, input A should be a s x s matrix;
% b and c should are column vectors of length s

%==========================================================================

tol = 1.e-15;
s = length(b); % number of stages
unity = ones(s,1);
p = 0;

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

% order 6
a(18) = b'*diag(c)^4*c;
a(19) = b'*diag(c)^3*A*c;
a(20) = b'*diag(c)*diag(A*c)*A*c;
a(21) = b'*diag(c)^2*A*diag(c)*c;
a(22) = b'*diag(A*c)*A*diag(c)*c;
a(23) = b'*diag(c)*A*diag(c)^2*c;
a(24) = b'*A*diag(c)^3*c;
a(25) = b'*diag(c)^2*A^2*c;
a(26) = b'*diag(A*c)*A^2*c;
a(27) = b'*diag(c)*A*diag(c)*A*c;
a(28) = b'*A*diag(c)^2*A*c;
a(29) = b'*A*diag(A*c)*A*c;
a(30) = b'*diag(c)*A^2*diag(c)*c;
a(31) = b'*A*diag(c)*A*diag(c)*c;
a(32) = b'*A^2*diag(c)^2*c;
a(33) = b'*diag(c)*A^3*c;
a(34) = b'*A*diag(c)*A^2*c;
a(35) = b'*A^2*diag(c)*A*c;
a(36) = b'*A^3*diag(c)*c;
a(37) = b'*A^4*c;

%==========================================================================

%% Order conditions

% order 1
t(1) = a(1) - 1;
if max(abs(t))<tol
  p = 1;
end

% order 2
t = a(2) - 1/2;

if abs(t) < tol && p == 1
    p = 2;
end

% order 3
t(1) = a(3) - 1/3;
t(2) = a(4) - 1/6;

if max(abs(t)) < tol && p == 2
    p = 3;
end

% order 4
t(1) = a(5) - 1/4;
t(2) = a(6)- 1/8;
t(3) = a(7)- 1/12;
t(4) = a(8) - 1/24;

if max(abs(t)) < tol && p == 3
    p = 4;
end

% order 5
t(1) = a(9) - 1/5;
t(2) = a(10) - 1/10;
t(3) = a(11) - 1/15;
t(4) = a(12) - 1/30;
t(5) = a(13)- 1/20;
t(6) = a(14) - 1/20;
t(7) = a(15)- 1/40;
t(8) = a(16) - 1/60;
t(9) = a(17) - 1/120;

if max(abs(t)) < tol && p == 4
    p = 5;
end

% order 6
t(1) = a(18) - 1/6;
t(2) = a(19) - 1/12;
t(3) = a(20) - 1/24;
t(4) = a(21) - 1/18;
t(5) = a(22) - 1/36;
t(6) = a(23) - 1/24;
t(7) = a(24) - 1/30;
t(8) = a(25) - 1/36;
t(9) = a(26) - 1/72;
t(10) = a(27) - 1/48;
t(11) = a(28) - 1/60;
t(12) = a(29) -1/120;
t(13) = a(30) - 1/72;
t(14) = a(31) - 1/90;
t(15) = a(32) - 1/120;
t(16) = a(33) - 1/144;
t(17) = a(34) - 1/180;
t(18) = a(35) - 1/240;
t(19) = a(36) - 1/360;
t(20) = a(37) - 1/720;

if max(abs(t)) < tol && p == 5
    p = 6;
    disp('This method has order at least 6');
end

end
