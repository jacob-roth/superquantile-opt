function [solution_str,outargstruc_arr] = PSG_solver(C,c,A,b,r,k,lb,ub,SOLVER)
% function [solution_str, problem_statement, outargstruc_arr] = PSG_solver(C,c,A,b,r,k,lb,ub,SOLVER)
%Initialize:
L = length(A);
[m,n] = size(A{1});
tau = 1 - k/m;
alpha = tau; % not 1-tau

%Pack data to PSG structure:
if sum(C) > 0
    disp('QUADRATIC objective');
    if issparse(C)
        iargstruc_arr(1) = pmatrix_pack('pmatrix_obj_data',0.5*C,[],c,[]); % 0.5x'*C*x + c'*x; NEED to call it 'matrix_xyz'
    else
        iargstruc_arr(1) = matrix_pack('matrix_obj_data',0.5*C,[],c,[]); % 0.5x'*C*x + c'*x; NEED to call it 'matrix_xyz'
    end
else
    disp('LINEAR objective');
    iargstruc_arr(1) = matrix_pack('matrix_obj_data',c',[],0.0,[]); % c'*x; NEED to call it 'matrix_xyz'
end
iargstruc_arr(2) = point_pack('point_lowerbounds',lb); % x >= lb; NEED to call it 'point_lowerbounds' for Box
iargstruc_arr(3) = point_pack('point_upperbounds',ub); % x <= ub; NEED to call it 'point_upperbounds' for Box
for ell = 1:L
  iargstruc_arr(3+ell) = matrix_pack(interpolate_string('matrix_con_data{ell}'),-A{ell},[],b{ell},[]); % A^ell x + b^ell; ; NEED to call it 'matrix_xyz'
end
 
%Constraints:
% con_ell = @(ell)(interpolate_string('Constraint: <= r_k{ell}\n  cvar_risk(alpha{ell},matrix_con_data{ell})\n'));
con_ell = @(ell)(interpolate_string('Constraint: <= r_k{ell},linearize=0\n  cvar_risk(alpha{ell},matrix_con_data{ell})\n'));
cvar_cons = '';
for ell = 1:L
  cvar_cons = strcat(cvar_cons,con_ell(ell));
end
cvar_cons = cvar_cons(1:end-2); % trim '\n'

if sum(C) > 0
    if issparse(C)
        %Define problem statement (quadratic:sparse):
        problem_statement = sprintf('%s\n',...
        'minimize',...
        'Objective: quadratic = 1',...
        '  quadratic(pmatrix_obj_data)',...
        sprintf(cvar_cons),...
        'Box_of_Variables: lowerbounds = point_lowerbounds, upperbounds = point_upperbounds',...
        'Solver: SOLVER, precision = 8, stages = 30, timelimit = 3600, print = 1'...
        ); % default is 6 stages per refman (http://www.aorda.com/html/PSG_Help_HTML/index.html?solver_in_text_format.htm) but often solves to insufficient precision; doesn't affect CVaR group
    else
        %Define problem statement (quadratic:full):
        problem_statement = sprintf('%s\n',...
        'minimize',...
        'Objective: quadratic = 1',...
        '  quadratic(matrix_obj_data)',...
        sprintf(cvar_cons),...
        'Box_of_Variables: lowerbounds = point_lowerbounds, upperbounds = point_upperbounds',...
        'Solver: SOLVER, precision = 8, stages = 30, timelimit = 3600, print = 1'...
        ); % default is 6 stages per refman (http://www.aorda.com/html/PSG_Help_HTML/index.html?solver_in_text_format.htm) but often solves to insufficient precision
    end
else
    %Define problem statement (linear):
    problem_statement = sprintf('%s\n',...
    'minimize',...
    'Objective: linearize = 1',...
    '  linear(matrix_obj_data)',...
    sprintf(cvar_cons),...
    'Box_of_Variables: lowerbounds = point_lowerbounds, upperbounds = point_upperbounds',...
    'Solver: SOLVER, precision = 8, stages = 30, timelimit = 3600, print = 1'...
    ); % default is 6 stages per refman (http://www.aorda.com/html/PSG_Help_HTML/index.html?solver_in_text_format.htm) but often solves to insufficient precision
end

%Change parameters with their values in problem statement:
for ell = 1:L
  problem_statement = strrep(problem_statement,interpolate_string('r_k{ell}'),num2str(r(ell)/k(ell),'%0.16f'));
  problem_statement = strrep(problem_statement,interpolate_string('alpha{ell}'),num2str(alpha(ell),'%0.16f'));
end
problem_statement = strrep(problem_statement,'SOLVER',SOLVER); %Solvers: VAN,TANK,CAR,BULDOZER,VANGRB,CARGRB,HELI

disp(problem_statement);

%Optimize problem using mpsg_solver function:
[solution_str, outargstruc_arr] = mpsg_solver(problem_statement, iargstruc_arr);
results = psg_convert_outargstruc(solution_str,outargstruc_arr);
outargstruc_arr = results.data;
solution_str = results.string;