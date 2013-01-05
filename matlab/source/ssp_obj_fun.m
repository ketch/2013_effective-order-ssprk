function [r,grad] = ssp_obj_fun(x)
% Objective function and its gradient for the optimization problem
%
% Used by effective_ssp.m

%==========================================================================

r = x(end);
grad = zeros(size(x));
grad(end) = 1;

end