function [f,gradf] = processors_obj_fun(x)
% Objective function and its gradient for the optimization problem
%
% Used by pre_post_processors.m

%==========================================================================

f = x(end);
gradf = zeros(length(x),1);
gradf(end) = 1;

end