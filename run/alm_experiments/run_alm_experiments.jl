ENV["GRB_LICENSE_FILE"] = "C:\\Users\\roth0674\\gurobi.lic"
if Sys.iswindows()
  # cd("C:\\Users\\roth0674\\Documents\\GitHub\\cc-non\\code\\run\\alm_experiments\\")
  cd("Y:/git/cc-non/code/run/alm_experiments/")
end
import Pkg; Pkg.activate("."); Pkg.instantiate()
using Distributed
# addprocs(4)

@everywhere begin
import Pkg; Pkg.activate("."); Pkg.instantiate()
using JuMP, Random
include("../../src/SSNALM.jl")
include("helper_alm_experiments.jl")
if Sys.iswindows()
  # global const DATAPATH = "C:\\Users\\roth0674\\alm_results\\"
  # global const psg_problemdata = "C:/Users/roth0674/Documents/GitHub/cc-non/code/run/alm_experiments/psg_problemdata/"
  # global const psg_problemoutput = "C:/Users/roth0674/Documents/GitHub/cc-non/code/run/alm_experiments/psg_problemoutput/"
  global const DATAPATH = "Y:/alm_results/"
  # global const psg_problemdata = "Y:/git/cc-non/code/run/alm_experiments/psg_problemdata/"
  # global const psg_problemoutput = "Y:/git/cc-non/code/run/alm_experiments/psg_problemoutput/"
  global const psg_problemdata = "C:/Users/roth0674/mytmp/psg_problemdata/"
  global const psg_problemoutput = "C:/Users/roth0674/mytmp/psg_problemoutput/"
  
  using MATLAB
  MATLAB_session = MSession() # `close(MATLAB_session)` to close it
else
  global const DATAPATH = "/home/roth0674/drive/alm_results/"
end
GC.enable_logging(false)

# global constants
global const maxtime = 3600
global const tol_hi = 1e-8
global const tol_lo1 = 1e-3
global const tol_lo2 = 1e-6
global const nthreads = 1

# alm parameters
lin_solver = (:M, 1)
AP = ALMParameters()
AP.maxait = 1_000
AP.stoptol = tol_hi
AP.stoptol_i1 = tol_lo1
AP.stoptol_i2 = tol_lo2
AP.maxn = 50_000
AP.maxtime = maxtime
BLAS.set_num_threads(nthreads)

# ssn parameters
SP = SSNParameters(lin_solver)
SP.maxnit = 200
SP.stagit = 10
SP.maxlit = 50
SP.tolsub = 1e-7
SP.ls_tau = 0.7
SP.ls_c1 = 1e-4
SP.ls_c2 = 0.9
SP.ls_method = :sw # {:bt, :sw, :ww}

# gurobi/osqp options
options = default_options()
options[:relgap_tol_hi] = tol_hi
options[:pinfeas_tol_hi] = tol_hi
options[:dinfeas_tol_hi] = tol_hi
options[:relgap_tol_lo] = tol_lo1
options[:pinfeas_tol_lo] = tol_lo1
options[:dinfeas_tol_lo] = tol_lo1
options[:maxtime] = maxtime
options[:nthreads] = nthreads
options[:presolve] = 0
end # everywhere

#
# warmup
#

writeout = false
for objtype in [1,2]
  pcts = [0.01]
  Ls = [1]
  nrep = 1
  seed = 1234567
  mtd = [:alm]#, :grb_oa, :grb, :psg, :osqp]

  # dimensions
  base_m = 10
  lo_m = 3
  hi_ms = [3]
  base_n = 10
  lo_n = 2
  hi_ns = [2]
  mn_step = 1
  #=
  # base_m = 10
  # lo_m = 5-2
  # hi_ms = [6-2,  6-2]
  # base_n = 10
  # lo_n = 2
  # hi_ns = [3,  3]
  # mn_step = 1
  =#
  outnames = []
  for pct in pcts
    for (i, L) in enumerate(Ls)
      mns = antidiagonal(
        base_m, lo_m, hi_ms[i], mn_step,
        base_n, lo_n, hi_ns[i], mn_step,
      )
      scale_row = false
      scale_col = true
      sigma0 = 1.0
      tau0 = 1.0
      t = 0.0 # difficulty: 1 = easy and 0 = hard
      TI0 = TestInstance(seed, objtype, 0, [0], [pct for l in 1:L], scale_col, scale_row, t, sigma0, tau0, AP, SP)
      runtime = @elapsed begin outname = test_alm_timing(
        pct, L, mtd, DATAPATH, mns[1], nrep,
        TI0,
        options,
        Float64, writeout
      )
      push!(outnames, outname)
      end # time
    end
  end
end

#
# experiments
#

for objtype in [2, 1]
  pcts = [0.01, 0.1]
  Ls = [1, 10]
  nrep = 1
  seed = 1234567
  # mtd = [:alm, :grb_oa, :grb, :psg, :osqp]
  # mtd = [:alm, :psg, :grb, :grb_oa]
  # mtd = [:alm, :grb]
  mtd = [:alm]

  # dimensions
  # base_m = 10
  # lo_m = [5, 4].-1
  # hi_ms = [6,  5].-1
  # base_n = 10
  # lo_n = [2, 2].-1
  # hi_ns = [3,  3].-1
  # mn_step = 1
  base_m = 2
  lo_ms = [14, 11]
  hi_ms = [20,  17]
  base_n = 2
  lo_ns = [7, 7]
  hi_ns = [13,  13]
  mn_step = 1
  outnames = []
  for pct in pcts
    for (i, L) in enumerate(Ls)
      mns = antidiagonal(
        base_m, lo_ms[i], hi_ms[i], mn_step,
        base_n, lo_ns[i], hi_ns[i], mn_step,
      )
      scale_row = false
      scale_col = true
      sigma0 = 1.0
      tau0 = 1.0
      t = 0.0 # difficulty: 1 = easy and 0 = hard
      TI0 = TestInstance(seed, objtype, 0, [0], [pct for l in 1:L], scale_col, scale_row, t, sigma0, tau0, AP, SP)
      runtime = @elapsed begin outname = test_alm_timing(
        pct, L, mtd, DATAPATH, mns[1], nrep,
        TI0,
        options,
        Float64
      )
      push!(outnames, outname)
      end
    end
  end
end