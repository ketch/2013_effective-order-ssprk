function obs_SSP_coef = table_obs_SSP_coef(problem,initial_data,x_left, ...
    x_right,t_final,gridpts)
% Creates table of observed SSP coefficiets for various problems
%
% Input variale meanings:
%    problem      --- 'advection','burgers','burgers2','buckley_leverett'
%    initial_data --- 'cont','discont', continuous and discontinuous data
%    x_left       --- left spatial point
%    x_right      --- right spatial point
%    t_final      --- final time
%    gridpts      --- number of spatial grid points
%
% Output variables:
%   obs_SSP_coef  --- table with obeserved SSP coefficients

%==========================================================================

format long;
clc; close all;

%==========================================================================

%% Matrix for observed SSP coefficients 
obs_SSP_coef = cell(7,13);
obs_SSP_coef{1,1} = 'q'; obs_SSP_coef{1,2} = 'p';
for i = 3:13
    obs_SSP_coef{1,i} = i-2;
end
obs_SSP_coef{2,1} = 2; obs_SSP_coef{2,2} = 2;
obs_SSP_coef{3,1} = 3; obs_SSP_coef{3,2} = 2;
obs_SSP_coef{4,1} = 3; obs_SSP_coef{4,2} = 3;
obs_SSP_coef{5,1} = 4; obs_SSP_coef{5,2} = 2;
obs_SSP_coef{6,1} = 4; obs_SSP_coef{6,2} = 3;
obs_SSP_coef{7,1} = 4; obs_SSP_coef{7,2} = 4;

SSP_coef = obs_SSP_coef;

i = 2;
for effective_order = 2:4
    for order = 2:effective_order
        for stages = effective_order+(order==4):11
            [sigma,r] = observed_ssp_coef(stages,effective_order,order, ...
                problem,initial_data,x_left,x_right,t_final,gridpts);
            obs_SSP_coef{i,stages+2} = sprintf('%1.2f(%2.f%%)',sigma, ... 
                sigma/r*100-100);
            SSP_coef{i,stages+2} = sprintf('%1.2f',r);
        end
    i = i+1;    
    end
end

end