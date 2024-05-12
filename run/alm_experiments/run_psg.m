% setup
addpath('C:\Aorda\PSG\MATLAB\Tools');
addpath('C:\Aorda\PSG\MATLAB');
addpath('C:\Aorda\PSG\lib');
% addpath('C:\Users\roth0674\Documents\GitHub\cc-non\code\src');
% addpath('C:\Users\roth0674\Documents\GitHub\cc-non\code\run\alm_experiments\');
% datapath = 'C:\Users\roth0674\Documents\GitHub\cc-non\code\run\alm_experiments\psg_problemdata\';
% outputpath = 'C:\Users\roth0674\Documents\GitHub\cc-non\code\run\alm_experiments\psg_problemoutput\';
% datapath = 'Y:\git\cc-non\code\run\alm_experiments\psg_problemdata\'; % too small storage
% outputpath = 'Y:\git\cc-non\code\run\alm_experiments\psg_problemoutput\'; % too small storage
addpath('Y:\git\cc-non\code\src');
addpath('Y:\git\cc-non\code\run\alm_experiments\');
addpath('C:\Users\roth0674\mytmp\');
datapath = 'C:\Users\roth0674\mytmp\psg_problemdata\';
outputpath = 'C:\Users\roth0674\mytmp\psg_problemoutput\';

% print
disp('initialized paths; calling PSG_solver_read');

% solve
% [solution_str,problem_statement,outargstruc_arr] = PSG_solver_read_sparse(datapath);
[solution_str,outargstruc_arr] = PSG_solver_read_sparse(datapath);

% print
disp('called solver; processing output');

% write
output_structure = tbpsg_solution_struct(solution_str, outargstruc_arr);
writematrix(output_structure.objective, strcat(outputpath,'pobj.csv'));
writematrix(output_structure.gap, strcat(outputpath,'gap.csv'));
writematrix(output_structure.time(1), strcat(outputpath,'load_time.csv'));
writematrix(output_structure.time(2), strcat(outputpath,'preprocess_time.csv'));
writematrix(output_structure.time(3), strcat(outputpath,'solve_time.csv'));
writematrix(output_structure.status{1}, strcat(outputpath,'status.csv'));
writematrix(output_structure.point_constraint_data,strcat(outputpath,'mks.csv'));
point_data = tbpsg_optimal_point_data(solution_str, outargstruc_arr);
writematrix(output_structure.optimal_point_data ,strcat(outputpath,'xbar.csv')); % just get directly from output_structure
% disp(point_data);
disp(strcat(outputpath,'xbar.csv'));
writematrix(point_data,strcat(outputpath,'xbar.csv'));
% exit