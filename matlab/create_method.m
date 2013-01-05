function create_method(stages,effective_order,order)
% Creates and stores ESSPRK methods is a .mat file.
%
% Input variable meanings:
%   stages          --- method's stages
%   effective_order --- method's effective order
%   order           --- method's classical order 

%==========================================================================

% Set path
cd Source

% create method
[R,Rb,Rc,M,b,c,T,Tb,Tc,r1,r,r2,effr1,ceff,effr2,P,Pb,Pc,order_P]...
            = get_method(stages,effective_order,order);
M,b,c
pause;
% check order of P = TMR
if order_P ~= effective_order
    sprintf(['RMT method with M = ESSPRK%d%d%d is not of desired ' ...
            'order %d'],stages, effective_order, order, effective_order)
    
    % store scheme
    file = sprintf('../Methods/ESSPRK%d%d%d.mat',stages, ... 
        effective_order,order);
    save(file,'R','Rb','Rc','M','b','c','T','Tb','Tc','r1','r','r2')
end

% Set path
cd ..

end