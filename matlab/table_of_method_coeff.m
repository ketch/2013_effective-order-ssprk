function table_of_method_coeff(stages,effective_order,order)
% Creates a three-column matrix with the EESSPRK scheme coefficients
% for practical usage in presentations.
% First column are coeff of starting method R, second column are 
% coefficients of main method M and last colum those of stopping method T.
%
% Input variable meanings:
%   stages          --- method's stages
%   effective_order --- method's effective order
%   order           --- method's classical order 

%==========================================================================

% Get coefficients
method = sprintf('Methods/ESSPRK%d%d%d.mat',stages,effective_order,order);
load(method);

% Coefficients of starting method
Rv = reshape(transpose(R),length(R)^2,1);
Rv(Rv==0) = [];
Rv = [Rv; Rb];

% Coefficients of main method
Mv = reshape(transpose(M),length(M)^2,1);
Mv(Mv==0) = [];
Mv = [Mv; b];

% Coefficients of stopping method
Tv = reshape(transpose(T),length(T)^2,1);
Tv(Tv==0) = [];
Tv = [Tv; Tb];

% Resulting matrix
v = cell(length(Rv),1);
j = 2;
k = 1;
for i = 1:length(Rv)
    if i <= length(Mv)
        if mod(k,j) == 0
            k = 1;
            j = j + 1;
        end
        if j <= length(b)
            v{i} = sprintf(['$a_{%d,%d}=%1.15f$ & $a_{%d,%d}=%1.15f$ ' ...
            '& $a_{%d,%d}=%1.15f$ \\\\'],j,mod(k,j),Rv(i),j,mod(k,j), ...
            Mv(i),j,mod(k,j),Tv(i));
            disp(v{i});
        else
            v{i} = sprintf(['$a_{%d,%d}=%1.15f$ & $b_%d=%1.15f$ & ' ...
                '$b_%d=%1.15f$ \\\\'],j,mod(k,j),Rv(i),mod(k,j),Mv(i), ...
                mod(k,j),Tv(i));
            disp(v{i});
        end
        k = k + 1;
    else
        v{i} = sprintf('$b_%d=%1.15f$ & $ $ & $ $ \\\\',i-length(Mv), ...
            Rv(i));
        disp(v{i});
    end
end
    
end
