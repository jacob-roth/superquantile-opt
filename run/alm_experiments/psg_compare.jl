import Pkg; Pkg.activate("."); Pkg.instantiate()
using JuMP, Random
include("../../src/SSNALM.jl")
include("helper_alm_experiments.jl")
global const DATAPATH = "/home/roth0674/drive/alm_results/"
GC.enable_logging(false)
BLAS.get_num_threads()
BLAS.set_num_threads(8)
Tf = Float64

#
# simple problem
#

# alm parameters
lin_solver = (:M, 1)
AP = ALMParameters()
AP.maxait = 1_000
AP.stoptol = 1e-6

# ssn parameters
SP = SSNParameters(lin_solver)
SP.maxnit = 200
SP.maxlit = 50
SP.tolsub = 1e-7
SP.ls_tau = 0.7
SP.ls_c1 = 1e-4
SP.ls_c2 = 0.9
SP.ls_method = :sw # {:bt, :sw, :ww}

# gurobi/osqp options
options = default_options()
options[:nthreads] = BLAS.get_num_threads()

# instance
pct = 0.01
L = 2
objtype = 2
TI = TestInstance(12345, objtype, 0, [0], [pct for l in 1:L], false, false, 0.0, -1.0, 0.0, AP, SP)
m = 3
n = 3
TI.m = [m for l in 1:L]
TI.n = n

# generate problem
Random.seed!(1234567)
# A = [Matrix(1.0I(3)), Matrix(1.0I(3))]
A = [randn(3,3), randn(3,3)]
b = [-2 .* [2.;-0.1;1.], -2 .* [2.;-0.1;1.]]
k = [1, 1]
m = [3, 3]
r = [-0.75; -0.75]
C = 1.0I(3)
Cinv = inv(C)
c = -1 .* [2.; 3.; 4.]
lb = [-1.; -1.; -1.]
ub = [1.; 1.; 1.]
lbidx = [true; true; true]
ubidx = [true; true; true]
Bt = Function[]
for l in eachindex(m)
  Btl = Atmatinv_action(m[l], k[l])
  push!(Bt, Btl)
end
CI = ConstraintInstance(
  TI.seed,
  TI.n, TI.m, L,
  A, b, k, r,
  Bt, lb, ub, lbidx, ubidx, false
)
OI = ObjectiveInstance(
  TI.seed, TI.objtype,
  C, Cinv,
  c
)
PI = ProblemInstance(CI, OI)
obj = FunctionValue(Tf, TI.n, TI.objtype) # 2==diagonal
if !isnothing(PI.C)
  obj.H .= PI.C
else
  obj.g .= PI.c
end

# solve gurobi
solve_grbps = @elapsed res = solve_cvar_gurobi(PI, options)
calcs_manual!(
  res, obj, PI,
  res[:x], res[:y], res[:z], res[:t],
  res[:lambda], res[:mu]
)
res