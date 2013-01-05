function save_method_to_file(stages,effective_order,order)
% Loads method from a .mat file and saves it to simple text file
%
% Input variable meanings:
%   stages          --- method's stages
%   effective_order --- method's effective order
%   order           --- method's classical order
%

%==========================================================================

%% Set method
method = sprintf('Methods/ESSPRK%d%d%d.mat',stages,effective_order, ... 
    order);
load(method)


%==========================================================================

%% Write to .txt file

% open file
fileID = fopen('method.txt','w');

% write starting method
coef = [R(:); Rb; r1];
s = length(Rb);
fprintf(fileID,'%1s\n','R');
fprintf(fileID,'%16.15f ',coef(1:s^2));
fprintf(fileID,'\n\n%1s\n','Rb');
fprintf(fileID,'%16.15f ',coef(s^2+1:s*(s+1)));
fprintf(fileID,'\n\n%1s\n','r1');
fprintf(fileID,'%16.15f\n\n\n',coef(end));

% write starting method
coef = [M(:); b; r];
s = length(b);
fprintf(fileID,'%1s\n','M');
fprintf(fileID,'%16.15f ',coef(1:s^2));
fprintf(fileID,'\n\n%1s\n','b');
fprintf(fileID,'%16.15f ',coef(s^2+1:s*(s+1)));
fprintf(fileID,'\n\n%1s\n','r');
fprintf(fileID,'%16.15f\n\n\n',coef(end));

% write starting method
coef = [T(:); Tb; r2];
s = length(Tb);
fprintf(fileID,'%1s\n','T');
fprintf(fileID,'%16.15f ',coef(1:s^2));
fprintf(fileID,'\n\n%1s\n','Tb');
fprintf(fileID,'%16.15f ',coef(s^2+1:s*(s+1)));
fprintf(fileID,'\n\n%1s\n','r2');
fprintf(fileID,'%16.15f',coef(end));

% close file
fclose(fileID);

end