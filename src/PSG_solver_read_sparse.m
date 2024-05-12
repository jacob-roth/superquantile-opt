function [solution_str,outargstruc_arr] = PSG_solver_read_sparse(datapath)
% function [solution_str,problem_statement,outargstruc_arr] = PSG_solver_read_sparse(datapath)

% Setup
L = readmatrix(strcat(datapath,'L.csv'));
SOLVER = fscanf(fopen(strcat(datapath,'SOLVER.csv')),'%s');
fclose('all');

% Objective
c = readmatrix(strcat(datapath,'cl.csv'));
n = length(c);
C = spalloc(n,n,n);
Cdiag = readmatrix(strcat(datapath,'Cq_diag.csv'));
for i = 1:n
    C(i,i) = Cdiag(i);
end
disp('populated Cdiag');

% Constraints
A = {};
b = {};
for ell = 1:L
    A{ell} = readmatrix(strcat(datapath,'A',num2str(ell),'.csv'));
    b{ell} = readmatrix(strcat(datapath,'b',num2str(ell),'.csv'));
end
r = readmatrix(strcat(datapath,'r.csv'));
k = readmatrix(strcat(datapath,'k.csv'));
lb = readmatrix(strcat(datapath,'lb.csv'));
ub = readmatrix(strcat(datapath,'ub.csv'));

%----------------- start optional: write to txt file -----------------
if 0==1
    disp('SAVING txt FILES');
    % Header: write to txt files
    varnames = {};
    for i = 1:n
        varnames{i} = strcat('x',num2str(i));
    end
    matrixheader = varnames;
    matrixheader{n+1} = 'scenario_benchmark'
    pointheader = {'component_name','value'}
    if issparse(C)
        [i,j,val] = find(C);
        iargstruc_arr(1) = pmatrix_pack('pmatrix_obj_data',0.5*C,[],c,[]); % 0.5x'*C*x + c'*x; NEED to call it 'matrix_xyz'
        fid = fopen(strcat(datapath,iargstruc_arr(1).name,'.txt'),'wt');
        fprintf(fid,'%s\t',string(matrixheader));
        fprintf(fid,'\n');
        for ctr = 1:n
            fprintf(fid,'%d\t',int64(i(ctr)));
            fprintf(fid,'%d\t',int64(j(ctr)));
            fprintf(fid,'%f\t',0.5*val(ctr));
            fprintf(fid,'\n');
        end
        for ctr = 1:n
            fprintf(fid,'%d\t',int64(ctr));
            fprintf(fid,'%d\t',int64(n+1));
            fprintf(fid,'%f\t',c(ctr));
            fprintf(fid,'\n');
        end
        fclose(fid);
    else
        iargstruc_arr(1) = matrix_pack('matrix_obj_data',0.5*C,[],c,[]); % 0.5x'*C*x + c'*x; NEED to call it 'matrix_xyz'
        writematrix([string(matrixheader); iargstruc_arr(1).data],strcat(datapath,iargstruc_arr(1).name,'.txt'),'Delimiter','tab')
    end
    iargstruc_arr(2) = point_pack('point_lowerbounds',lb); % x >= lb; NEED to call it 'point_lowerbounds' for Box
    writematrix([string(pointheader); [string(varnames)', iargstruc_arr(2).data']],strcat(datapath,iargstruc_arr(2).name,'.txt'),'Delimiter','tab')
    iargstruc_arr(3) = point_pack('point_upperbounds',ub); % x <= ub; NEED to call it 'point_upperbounds' for Box
    writematrix([string(pointheader); [string(varnames)', iargstruc_arr(3).data']],strcat(datapath,iargstruc_arr(3).name,'.txt'),'Delimiter','tab')
    for ell = 1:L
        iargstruc_arr(3+ell) = matrix_pack(strcat('matrix_con_data',num2str(ell)),-A{ell},[],b{ell},[]); % A^ell x + b^ell; ; NEED to call it 'matrix_xyz'
        writematrix([string(matrixheader); iargstruc_arr(3+ell).data],strcat(datapath,iargstruc_arr(3+ell).name,'.txt'),'Delimiter','tab')
    end
end
%----------------- end optional: write to txt file -----------------

display(strcat('calling PSG_solver with:',SOLVER));
[solution_str, outargstruc_arr] = PSG_solver(C,c,A,b,r,k,lb,ub,SOLVER);

%----------------- start optional: write to txt file -----------------
if 0==1
    fid = fopen(strcat(datapath,'problem_mks.txt'),'wt');
    fprintf(fid,'%s',problem_statement);
    fclose(fid);
end
%----------------- end optional: write to txt file -----------------