function [SSP_coef, log] = create_methods(effective_order,order)
% Create explicit ESSPRK methods and stores effective SSP coefficients
% to check optimality.
%
% Input variable meanings:
%   effective_order --- method's effective order
%   order           --- method's classical order 
%
% Output variables:
%   SSP_coef        --- table with effective SSP coefficients
%   log             --- log file of errors if a method is not of the 
%                       desired order

%==========================================================================

% Set path
cd Source

%==========================================================================

%% Matrix for SSP coefficients 
SSP_coef = cell(7,13);
SSP_coef{1,1} = 'q'; SSP_coef{1,2} = 'p';
for i = 3:13
    SSP_coef{1,i} = i-2;
end
SSP_coef{2,1} = 2; SSP_coef{2,2} = 2;
SSP_coef{3,1} = 3; SSP_coef{3,2} = 2;
SSP_coef{4,1} = 3; SSP_coef{4,2} = 3;
SSP_coef{5,1} = 4; SSP_coef{5,2} = 2;
SSP_coef{6,1} = 4; SSP_coef{6,2} = 3;
SSP_coef{7,1} = 4; SSP_coef{7,2} = 4;
log = [];

%==========================================================================

%% Create effective order SSPRK methods

% ESSPRK(s,2,2)
if effective_order == 2 && order == 2
    for stages = 2:11
        
        % known SSP coefficients
        coef = [0.5 0.67 0.75 0.8 0.83 0.86 0.88 0.89 0.9 0.91];
        
        % get method
        flag = 1;
        while flag
            [R,Rb,Rc,M,b,c,T,Tb,Tc,r1,r,r2,effr1,ceff,effr2,P,Pb,Pc, ...
                order_P] = get_method(stages,effective_order,order);
            if abs(coef(stages-1) - ceff) < 1e-2
                flag = 0;
            end
        end
        
        % store scheme
        file = sprintf('../Methods/ESSPRK%d%d%d.mat',stages, ...
            effective_order,order);
        save(file,'R','Rb','Rc','M','b','c','T','Tb','Tc','r1','r','r2')
        
        % store SSP coefficient
        SSP_coef{2,stages+2} = sprintf('%1.2f',r);
        
        % check order of P = TMR
        if order_P ~= 2
            log = [log; sprintf(['RMT method with M = ESSPRK%d%d%d is ' ...
                'not of desired order %d'],stages, effective_order, ... 
                order, effective_order)];
        end
    end
end

% ESSPRK(s,3,2)
if effective_order == 3 && order == 2
    for stages = 3:11
        
        % known SSP coefficients
        coef = [0.33 0.50 0.53 0.59 0.61 0.64 0.67 0.68 0.69];
        
        % get method
        flag = 1;
        while flag
            [R,Rb,Rc,M,b,c,T,Tb,Tc,r1,r,r2,effr1,ceff,effr2,P,Pb,Pc, ...
                order_P] = get_method(stages,effective_order,order);
            if abs(coef(stages-2) - ceff) < 1e-2
                flag = 0;
            end
        end
        
        % store scheme
        file = sprintf('../Methods/ESSPRK%d%d%d.mat',stages, ... 
            effective_order,order);
        save(file,'R','Rb','Rc','M','b','c','T','Tb','Tc','r1','r','r2')
        
        % store SSP coefficient
        SSP_coef{3,stages+2} = sprintf('%1.2f',r);
        
        % check order of P = TMR
        if order_P ~= 3
            log = [log; sprintf(['RMT method with M = ESSPRK%d%d%d is ' ...
                'not of desired order %d'],stages, effective_order, ... 
                order, effective_order)];
        end
    end
end

% ESSPRK(s,3,3)
if effective_order == 3 && order == 3
    for stages = 3:11
        
        % known SSP coefficients
        coef = [0.33 0.50 0.53 0.59 0.61 0.64 0.67 0.68 0.69];
        
        % get method
        flag = 1;
        while flag
            [R,Rb,Rc,M,b,c,T,Tb,Tc,r1,r,r2,effr1,ceff,effr2,P,Pb,Pc, ...
                order_P] = get_method(stages,effective_order,order);
            if abs(coef(stages-2) - ceff) < 1e-2
                flag = 0;
            end
        end
        
        % store scheme
        file = sprintf('../Methods/ESSPRK%d%d%d.mat',stages, ... 
            effective_order,order);
        save(file,'R','Rb','Rc','M','b','c','T','Tb','Tc','r1','r','r2')
        
        % store SSP coefficient
        SSP_coef{4,stages+2} = sprintf('%1.2f',r);
        
        % check order of P = TMR
        if order_P ~= 3
            log = [log; sprintf(['RMT method with M = ESSPRK%d%d%d is ' ...
                'not of desired order %d'],stages, effective_order, ... 
                order, effective_order)];
        end
    end
end

% ESSPRK(s,4,2)
if effective_order == 4 && order == 2
    for stages = 4:11
        
        % known SSP coefficients
        coef = [0.22 0.39 0.44 0.50 0.54 0.57 0.60 0.62];
        
        % get method
        flag = 1;
        while flag
            [R,Rb,Rc,M,b,c,T,Tb,Tc,r1,r,r2,effr1,ceff,effr2,P,Pb,Pc, ...
                order_P] = get_method(stages,effective_order,order);
            if abs(coef(stages-3) - ceff) < 1e-2
                flag = 0;
            end
        end
        
        % store scheme
        file = sprintf('../Methods/ESSPRK%d%d%d.mat',stages, ... 
            effective_order,order);
        save(file,'R','Rb','Rc','M','b','c','T','Tb','Tc','r1','r','r2')
        
        % store SSP coefficient
        SSP_coef{5,stages+2} = sprintf('%1.2f',r);
        
        % check order of P = TMR
        if order_P ~= 4
            log = [log; sprintf(['RMT method with M = ESSPRK%d%d%d is ' ...
                'not of desired order %d'],stages, effective_order, ... 
                order, effective_order)];
        end
    end
end

% ESSPRK(s,4,3)
if effective_order == 4 && order == 3
    for stages = 10:10
        
        % known SSP coefficients
        coef = [0.19 0.37 0.43 0.50 0.54 0.57 0.60 0.62];
        
        % get method
        flag = 1;
        while flag
            [R,Rb,Rc,M,b,c,T,Tb,Tc,r1,r,r2,effr1,ceff,effr2,P,Pb,Pc, ...
                order_P] = get_method(stages,effective_order,order);
            if abs(coef(stages-3) - ceff) < 1e-2
                flag = 0;
            end
        end
        
        % store scheme
        file = sprintf('../Methods/ESSPRK%d%d%d.mat',stages, ... 
            effective_order,order);
        save(file,'R','Rb','Rc','M','b','c','T','Tb','Tc','r1','r','r2')
        
        % store SSP coefficient
        SSP_coef{6,stages+2} = sprintf('%1.2f',r);
        
        % check order of P = TMR
        if order_P ~= 4
            log = [log; sprintf(['RMT method with M = ESSPRK%d%d%d is ' ...
                'not of desired order %d'],stages, effective_order, ... 
                order, effective_order)];
        end
    end
end

% ESSPRK(s,4,4)
if effective_order == 4 && order == 4
    for stages = 11:11
        
        % known SSP coefficients
        coef = [0.30 0.38 0.47 0.52 0.54 0.60 0.59];
        
        % get method
        flag = 1;
        while flag
            [R,Rb,Rc,M,b,c,T,Tb,Tc,r1,r,r2,effr1,ceff,effr2,P,Pb,Pc, ...
                order_P] = get_method(stages,effective_order,order);
            if abs(coef(stages-4) - ceff) < 1e-2
                flag = 0;
            end
        end
        
        % store scheme
        file = sprintf('../Methods/ESSPRK%d%d%d.mat',stages, ... 
            effective_order,order);
        save(file,'R','Rb','Rc','M','b','c','T','Tb','Tc','r1','r','r2')
        
        % store SSP coefficient
        SSP_coef{7,stages+2} = sprintf('%1.2f',r);
        
        % check order of P = TMR
        if order_P ~= 4
            log = [log; sprintf(['RMT method with M = ESSPRK%d%d%d is ' ...
                'not of desired order %d'],stages, effective_order, ... 
                order, effective_order)];
        end
    end
end

% Set path
cd ..

end