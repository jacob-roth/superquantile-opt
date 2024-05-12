addpath('C:\Aorda\PSG\MATLAB\Tools');
addpath('C:\Aorda\PSG\MATLAB');
addpath('C:\Aorda\PSG\lib');
     

% A1 = [1 0 0; 0 1 0; 0 0 1];
% A2 = [1 0 0; 0 1 0; 0 0 1];
A1 = readmatrix('C:\Users\roth0674\A1.csv');
A2 = readmatrix('C:\Users\roth0674\A2.csv');
b1 = -2.*[2;-0.1;1];
b2 = -2.*[2;-0.1;1];
A = {A1, A2};
b = {b1, b2};
k = [1;1];
r = [-0.75;-0.75];
C = [1 0 0; 0 1 0; 0 0 1];
c = -1.*[2;3;4];
lb = [-1;-1;-1];
ub = [1;1;1];
SOLVER = 'VAN';
[solution_str,outargstruc_arr,xbar,ps,res] = PSG_solver(C,c,A,b,r,k,lb,ub,SOLVER);
xbar
