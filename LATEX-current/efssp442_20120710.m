cM = [0
0.730429885783319
0.644964638145795
1.000000000000000];

AM = [0 0 0 0;
0.730429885783319 0 0 0;
 0.251830917810810 0.393133720334985 0 0;
 0.141062771617064 0.220213358584678 0.638723869798257 0];

bM = [0.384422161080494 0.261154113377550 0.127250689937518 0.227173035604438]';

cR = [0
      0.545722177514735
0.842931687441527
0.574760809487828
0.980872743236632];

AR = [0 0 0 0 0;
      0.545722177514735 0 0 0 0;
0.366499989048164  0.476431698393363  0 0 0;
0.135697968350722 0.176400587890242 0.262662253246864 0 0;
0.103648417776838 0.134737771331049 0.200625899485633 0.541860654643112 0];

bR = [0.233699169638954 0.294263351266422 0.065226988215286 0.176168374199685 0.230642116679654]';

cT = [0
      0.509877496215340
      0.435774135529007
      0.933203341300203];

AT = [0 0 0 0;
      0.509877496215340 0 0 0;
      0.182230305923759 0.253543829605247 0 0;
      0.148498121305090 0.206610981494095 0.578094238501017 0];


bT = [0.307865440399752 0.171863794704750 0.233603236964822 0.286667527930676]';



% That's one big method for mankind
c = [cR; cM+1; cT+2];
b = [bR; bM; bT];

e4 = ones(4,1);
Z44 = zeros(4,4);
Z54 = zeros(5,4);

A = [ ...
    AR         Z54        Z54; ...
    e4 * bR'   AM         Z44; ...
    e4 * bR'   e4 * bM'   AT]

% one smaller step by a third
c = c / 3;
b = b / 3;
A = A / 3;


%% test M instead, 2nd-order as promised :-P
A = AM;  c = cM;  b = bM;
%A = AR;  c = cR;  b = bR;
%A = AT;  c = cT;  b = bT;

%% test RT
if (1==0)
c = [cR; cT+1];
b = [bR; bT];

e4 = ones(4,1);
Z44 = zeros(4,4);
Z54 = zeros(5,4);

A = [ ...
    AR         Z54; ...
    e4 * bR'   AT]

c = c / 2;
b = b / 2;
A = A / 2;
end


C = diag(c);
e = ones(length(c),1);

disp('Ae = c:')
max(abs(A * e - c))
disp('OC errors:')
% order 1
oc(1) = b' * e;
% order 2
oc(2) = b' * c;
% order 3
oc(3) = b' * c.^2;
oc(4) = b' * A * c;
% order 4
oc(5) = b' * c.^3;
oc(6) = b' * C * A * c;
oc(7) = b' * A * c.^2;
oc(8) = b' * A * A * c;

g = [1 1/2 1/3 1/6 1/4 1/8 1/12 1/24];

oc - g

if (1==0)
%beta2 = 0.0861063304588688/2

alpha = oc;
beta = zeros(size(alpha));
beta(2) = -1/6 + 1/2*alpha(3);
beta(3) = 1/12 - 1/2*alpha(3) + 1/3*alpha(5);
beta(4) = -1/24-1/3*alpha(5) + alpha(6);

% efOC 42
1/4 - alpha(3) + alpha(5) - 2*alpha(6) + alpha(7)
end


alpha(1) - rho(1)
alpha(2) + beta(2) - rho(2)
alpha(3) + beta(3) - rho(3)