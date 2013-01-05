function y_new = onestepRK(t,y_old,dt,A,b,c,problem,aux)
% Calculates the numerical solution over a single time-step using an
% RK method
%
% Used by test_method.m, solve_HCL.m, observed_ssp_coef.m,
% compute_dtFE_exp.m

%==========================================================================

s = length(b);
% Preallocation of step estimations
Y = zeros(s,length(y_old));

y_new = y_old;
for i=1:s
    Y(i,:) = y_old;
    for j=1:i-1
        Y(i,:) = Y(i,:) + dt.*A(i,j).*myfun(t+c(j)*dt,Y(j,:),problem,aux);
    end
    y_new = y_new + dt.*b(i).*myfun(t+c(i)*dt,Y(i,:),problem,aux);
end

end
