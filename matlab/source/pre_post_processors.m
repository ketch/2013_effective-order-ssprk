function [R,Rb,Rc,r1,effr1,T,Tb,Tc,r2,effr2] = pre_post_processors(M,b, ...
    c,r,effp,p,class)
% Find optimal SSPRK starting and stopping methods relative to the main 
% ESSPRK(s,q,p).
% The code works for explicit and singly implicit ESSPRK methods and 
% starting and stopping methods are explicit or singly implicit depending 
% on the main method.
%
% Input variable meanings:
%   M, b, c               -- main ESSPRK method's coefficients
%   r                     -- radius of absolute monotonicity 
%   s                     -- number of stages
%   effp                  -- effective order of accuracy
%   p                     -- order of accuracy
%   class                 -- explicit or singly implicit
%
% Output variables:
%   R,T, Rb,Tb, Rc,Tc     -- Method's coefficients
%   R,T are s x s
%   Rb,Tb re s x 1
%   Rc,Tc are s x 1
%   r1,r2 (scalar)        -- radius of absolute monotonicity
%   ceff1,ceff2 (scalar)  -- effective SSP coefficients for R and T
%
% Stored in a single vector x as:
% x = [R Rb' T Tb' -r1 -r2 -min(r1,r2)]
% Matrices are stored row-by-row
%
% Used by get_method.m

%==========================================================================

% Editable options:
solveorderconditions = 1; % starting guess satisfying equality constraints
abscissalessthanone = 1; % additional inequality such that absicca c <= 1
minr = 0.0; % Keep looking until a method with at least this value is found
maxr = -50.0; % Upper bund for radius of absolute monotonicity (r.a.m.)

%==========================================================================

% Input variables
if nargin == 6
    explicit = 1;
elseif strcmp(class,'IM')
    explicit = 0; % implicit method
end

Ms = length(b);
if explicit == 1 && Ms > 5 && effp == 4 && p == 2 && ~mod(sqrt(Ms-1),1)
    abscissalessthanone = 0;
end

%==========================================================================

%% Case effp = p
% Then method S in S^{-1}MS is the identity, hence R = T = M

if effp == p
    R = M; T = M;
    Rb = b; Tb = b;
    Rc = c; Tc = c;
    r1 = r; r2 = r;
    effr1 = r/length(b); effr2 = r/length(b);
    return;
end

%% Optimization routine

% Set optimization parameters:
opts = optimset('MaxFunEvals',1000000,'TolCon',1.e-14,'TolFun',1.e-14,...
    'TolX',1.e-14,'GradObj','on','MaxIter',500,'Diagnostics','on',...
    'Display','iter','Algorithm','sqp','RelLineSrchBnd',0.1,...
    'RelLineSrchBndDuration',100000000);
opts2 = optimset('Algorithm','levenberg-marquardt');
par = [explicit,abscissalessthanone];

%==========================================================================

% Set the number of stages for methods R and T
Rs = length(b)+1;
Ts = length(b);
if (effp~=p)
    Rs = Rs + 1;
end

% Set number of conditions on beta
s = 2*(effp==3) + 4*(effp==4);

% Include slack variable for objective function: min(r1,r2)
s = s + 1;

% Set the number of unknowns    -- n
% Variables for matrix of coefficients (strictly lower triangular if main 
% method is explicit, otherwise lower traingular if main method is 
% implicit) and for weights vector.
% The radii of absolute monotonicity and their minimum are the 3 entries in
% the unknowns vector.
explicitn = (Rs+1)*(Rs/2)+1 + (Ts+1)*(Ts/2)+1 + s;
implicitn = (Rs+3)*(Rs/2)+1 + (Ts+3)*(Ts/2)+1 + s;
n = explicit*explicitn + (1-explicit)*implicitn;

% Set the upper and lower bounds on the unknowns    -- ub, lb
lb = -5+zeros(1,n); lb(end-2) = maxr; lb(end-1) = maxr; lb(end) = maxr;
ub = 5+zeros(1,n); ub(end-2) = minr; ub(end-1) = minr; ub(end) = minr;

%==========================================================================

r1 = 0.0; 
r2 = 0.0;
info = -2;
while (info == -2 || r1 < r || r2 < r || info == 0)

    %Set initial guess
    x(1:n) = 1/(2*s)*(rand(1,n)); % currently only random guesses
    x(end-2:end) = [-r -r -r];
    
    %======================================================================

    % Find a feasible (for the order conditions) point to start
    % Note: code is much faster using this option
    if solveorderconditions == 1
        x = fsolve(@(x) oc_processors(x,Rs,Ts,M,b,c,effp,par(1)),x,opts2);
    end
  
    %======================================================================
  
    % The optimization call
    [X,FVAL,info] = fmincon(@processors_obj_fun,x,[],[],[],[],lb,ub, ...
        @(x) processors_nonlinear_con(x,Rs,Ts,M,b,c,effp,par),opts);
    r1 = -X(end-2); %Radius of absolute monotonicity of R
    r2 = -X(end-1); %Radius of absolute monotonicity of T
end 

%==========================================================================

%% Output

% Extract the Butcher coefficients from the solution X
if explicit == 1
    X1 = X(1:(Rs+1)*Rs/2);
    X2 = X(length(X1)+1:length(X1)+(Ts+1)*Ts/2);
else
    X1 = X(1:(Rs+3)*Rs/2);
    X2 = X(length(X1)+1:length(X1)+(Ts+3)*Ts/2);
end
[R,Rb] = unpack_RK(X1,Rs,explicit,0);
Rc = sum(R,2);
[T,Tb] = unpack_RK(X2,Ts,explicit,0);
Tc = sum(T,2);

% Effective SSP coefficients
effr1 = r1/Rs;
effr2 = r2/Ts;

end
