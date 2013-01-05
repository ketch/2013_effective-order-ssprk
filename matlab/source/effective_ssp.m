function [A,b,c,r,ceff,alpha,beta] = effective_ssp(s,effp,p,class)
% Find optimal effective order q SSPRK methods (ESSPRK(s,q,p)) with 
% s stages and classical order p.
% The code works for explicit and singly implicit methods.
%
% Input variable meanings:
%   s               -- number of stages
%   effp            -- effective order of accuracy
%   p               -- order of accuracy
%   class           -- explicit or singly implicit
%
% Output variables:
%   A, b, c         -- Method's coefficients
%   A is s x s
%   b is s x 1
%   c is s x 1
%   r (scalar)      -- radius of absolute monotonicity
%   ceff (scalar)   -- effective SSP coefficient (r/s)
%
% Stored in a single vector x as:
% x = [A b' -r]
% Matrices are stored row-by-row
%
% Special explicit cases: ESSPRK(v,2,2), ESSPRK(v^2,3,3), ESSPRK(v^2+1,4,2) 
%                         for which method's coefficients are known.
%
% Method is y = u^n + Dt*SUM(A*F(t_n + cDt,y^j))
%           u^n+1= u^n + Dt*SUM(b*F(t_n + cDt,y)
%
% Used by get_method.m

% TODO: Check order conditions of implicit methods for orders more than 5

%==========================================================================

% Editable options:
solveorderconditions = 1; % starting guess satisfying equality constraints
optspecialcase = 0; % use optimization for ESSPRK(n^2+1,4,2) case
shuoshercoefficients = 0; % opt variables set to Shu-Osher, default:Butcher
abscissalessthanone = 1; % additional inequality such that absicca c <= 1
minr = 0.0; % Keep looking until a method with at least this value is found
maxr = -50.0; % Upper bund for radius of absolute monotonicity (r.a.m.)

%==========================================================================

% Input variables
if nargin == 3
    explicit = 1;
elseif strcmp(class,'IM')
    explicit = 0; % implicit method
end

%==========================================================================

% Input checks
if p == 1 && effp > 1
    error(['Methods with effective order greater than one, cannot ', ... 
        'have classical order equal to one.'])
end
if explicit == 1 && (effp < p || s < p)
    error('Effective order must be at least equal to classical order.')
end
if explicit == 1 && s < effp
    error('Number of stages must be at least equal to effecrtive order.')
end

%==========================================================================

% Known coefficients cases
if explicit == 1 && optspecialcase == 0
    if s > 1 && effp == 2 && p == 2
        % load ESSPRK(v,2,2), v > 1
        [alpha,beta,A,b,c,r,ceff] = load_method(s,'ESSPRK2');
        return;
    elseif s > 4 && effp == 3 && p == 3 && ~mod(sqrt(s),1)
        % load ESSPRK(v^2,3,3), v > 2
        [alpha,beta,A,b,c,r,ceff] = load_method(s,'ESSPRK3');
        return;
    elseif s > 5 && effp == 4 && p == 2 && ~mod(sqrt(s-1),1)
        % load ESSPRK(v^2+1,4,2), v > 2
        [alpha,beta,A,b,c,r,ceff] = load_method(s,'ESSPRK4');
        return;
    end
end

%==========================================================================

%% Optimization routine

% Set optimization parameters:
opts = optimset('MaxFunEvals',1000000,'TolCon',1.e-14,'TolFun',1.e-14,...
    'TolX',1.e-14,'GradObj','on','MaxIter',500,'Diagnostics','on',...
    'Display','iter','Algorithm','sqp','RelLineSrchBnd',0.1,...
    'RelLineSrchBndDuration',100000000);
opts2 = optimset('Algorithm','levenberg-marquardt');
par = [explicit,shuoshercoefficients,abscissalessthanone];

%==========================================================================

% Set the number of unknowns    -- n
% Variables for matrix of coefficients (strictly lower triangular if method
% is explicit, otherwise lower traingular if method is implicit) and for 
% weights vector.
% The radius of absolute monotonicity is the last entry in the unknowns 
% vector.
n = explicit*(s+1)*(s/2) + (1-explicit)*(s+3)*(s/2) + 1;

% Set the upper and lower bounds on the unknowns    -- ub, lb
lb = -5.0 + zeros(1,n); lb(end) = maxr;
ub = 5.0 + zeros(1,n); ub(end) = minr;

%==========================================================================

r = -1.0;
info = -2;
while info == -2 || r < minr || info == 0
    
  % Set initial guess
  x(1:n) = 1/(2*s)*(rand(1,n)); % currently only random guesses
  if explicit == 1
      x(end) = -(1-solveorderconditions)*0.1 - ...
          solveorderconditions*(s-effp+1); % max r.a.m. (linear problems) 
  else
      x(end) = -0.1;
  end
  %========================================================================

  % Find a feasible (for the order conditions) point to start
  % Note: code is much faster using this option
  if solveorderconditions == 1
      x = fsolve(@(x) oc_effective_ssp(x,s,p,effp,par(1),par(2)),x,opts2);
  end
  
  %========================================================================
  
  % The optimization call
  [X,FVAL,info] = fmincon(@ssp_obj_fun,x,[],[],[],[],lb,ub,@(x) ...
      ssp_nonlinear_con(x,s,p,effp,par),opts);
  r = -FVAL;
  
end

%==========================================================================

%% Output

if shuoshercoefficients == 0
    % Extract the Butcher coefficients from the solution X
    [A,b] = unpack_RK(X,s,explicit,shuoshercoefficients);
    c = sum(A,2);
    
    % Construct Shu-Osher matrices
    [alpha,beta] = butcher2shuosher(A,b,r);
else
    % Extract the Shu-0sher coefficients from the solution X
    [alpha,beta] = unpack_RK(X,s,explicit,shuoshercoefficients);
    
    % Construct Butcher tableau
    [A,b] = shuosher2butcher(alpha,beta);
    c = sum(A,2);
end

% Effective SSP coefficient
ceff = r/s;

end
