function tv = tvnorm(y)
% Calculates the TV norm of a vector
%
% Used by test_method.m, solve_HCL.m, observed_ssp_coef.m,
% compute_dtFE_exp.m

%==========================================================================

tv = sum(abs(y(2:end) - y(1:end-1)));
tv = tv + abs(y(1) - y(end));

end