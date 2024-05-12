import Pkg; Pkg.activate("."); Pkg.instantiate()
using JuMP, Random, Parquet2, DataFrames, StatsBase
include("../../src/SSNALM.jl")
include("helper_taxi_experiments.jl")
GC.enable_logging(false)

# global constants
global const maxtime = 3600
global const tol_hi = 1e-8
global const tol_lo1 = 1e-3
global const tol_lo2 = 1e-6
global const TOL = 1e-8
global const nthreads = 1
global const DATAPATH = "/home/roth0674/drive/taxidata/parquet_combined2/"
global const DATAFILENAME = "taxi_weather_data_nomissing.parquet"
global const RESULTPATH = "/home/roth0674/drive/taxi_results/"

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
SP.stagit = 2
SP.maxnit = 200
SP.maxlit = 50
SP.tolsub = 1e-7
SP.ls_tau = 0.7
SP.ls_c1 = 1e-4
SP.ls_c2 = 0.9
SP.ls_method = :sw # {:bt, :sw, :ww}

# grb options
options = default_options()
# options[:pinfeas_tol_hi] = tol_hi
# options[:dinfeas_tol_hi] = tol_hi
# options[:relgap_tol_hi] = tol_hi
options[:pinfeas_tol_hi] = TOL # == 1e-4
options[:dinfeas_tol_hi] = TOL # == 1e-4
options[:relgap_tol_hi] = TOL # == 1e-4
options[:maxtime] = maxtime
options[:nthreads] = nthreads
options[:presolve] = 0

#
# taxi problem
#

XX, yy, XXnames = get_regression_data(DATAPATH, DATAFILENAME)
allcols = [
  "trip_miles",
  "trip_time",
  "base_passenger_fare",
  "tolls",
  "driver_pay",
  "delay_1",
  "delay_2",
  "avg_speed",
  "UBER",
  "PUborough_EWR",
  "PUborough_QNS",
  "PUborough_BNX",
  "PUborough_STI",
  "PUborough_BKN",
  "PUborough_UNK",
  "DOborough_EWR",
  "DOborough_QNS",
  "DOborough_BNX",
  "DOborough_STI",
  "DOborough_BKN",
  "DOborough_UNK",
  "shared_match_flag_Y",
  "shared_request_flag_Y",
  "wav_match_flag_Y",
  "wav_request_flag_Y",
  "tod_early",
  "tod_morning",
  "tod_evening",
  "tod_night",
  "tod_late",
  "visibility",
  "dew_point",
  "wind_speed",
  "wind_gust",
  "temp",
  "feels_like",
  "humidity",
  "pressure",
  "rain1h",
  "rain3h",
  "weekend",
  "fare_per_mile",
  "fare_per_min",
  "driver_share",
  "pay_per_mile",
  "pay_per_min",
  "congestion_per_min",
  "congestion_per_mile",
  "cloudy",
  "congestion_Y",
]
cols = [
  ### main ###
  "trip_time",
  "avg_speed",
  "base_passenger_fare",
  "congestion_per_min",
  "congestion_Y", 
  "shared_match_flag_Y",
  "shared_request_flag_Y",
  "wav_match_flag_Y",
  "wav_request_flag_Y",
  "tod_early",
  "tod_morning",
  "tod_evening",
  "tod_night",
  "tod_late",
  "tolls",
  "temp",
  "delay_1",
  "delay_2",
  ### others ###
  # "trip_time",
  # "temp",
  # "delay_1",
  # "delay_2",
  # "PUborough_EWR", # only ~100 "yes"
  # "DOborough_EWR", # only ~100 "yes"
]
col2idx = Dict()
for n in XXnames
  col2idx[n] = findfirst(n.==XXnames)
end
ZZ = XX[:, [col2idx[n] for n in cols]]

m = length(yy)
A = [-ZZ -ones(m)];
b = yy;
c = -vec(sum(A, dims=1)/m);
c[abs.(c) .< 1e-9] .= 0.0;
m, n = size(A)

#
# cold-start timing
#

global const TOL = 1e-8
pcts = [0.01, 0.05, 0.1, 0.25, 0.5]
for pct in pcts
  outname = "k_$(Int(pct*1000))pct__L_1"
  r = 0.0
  k = Int(ceil(pct*m))
  PI = ProblemInstance(
    -1, -1, n, [m], 1, [A], [b], [k], [r],
    Function[Atmatinv_action(m, k)],
    fill(-Inf, n), fill(+Inf, n), fill(false, n), fill(false, n), true,
    1, nothing, nothing, c[1:n], Dict()
  )

  Tf = Float64
  obj = FunctionValue(Tf, PI.n, PI.objtype) # 2==diagonal
  obj.g .= PI.c
  BLAS.set_num_threads(nthreads)
  x0 = zeros(Tf, PI.n);
  x0[end] = 1.0;
  lambda0 = [zeros(m[l]) for l in eachindex(m)];
  mux0 = zeros(Tf, PI.n);
  debug = false
  verb = 2
  sigma0 = 1e-5 #! worse: sigma0 = 1e-8
  tau0 = 1e-6   #! worse: tau0 = 1e-9
  
  lambda0 = [zeros(PI.m[l]) for l in eachindex(PI.m)];
  mux0 = zeros(Tf, PI.n);
  debug = false
  verb = 1
  tau0 = 1.0
  kgrid = nothing
  fullit = sum(PI.m)/sum(PI.L) <= 2PI.n
  SP.maxnit = 200
  SP.stagit = sum(PI.m)/PI.L > 4PI.n ? 10 : 100
  AP.BXratio = sqrt(sum(PI.m)/sum(PI.L)/PI.n)
  sigma0 = 1 / AP.BXratio
  tau0 = sigma0 / AP.BXratio
  # solve alm
  begin
    # GC.gc(); solve_alm = @elapsed res, _, _ = alm(x0, lambda0, mux0, sigma0, tau0, PI, AP, SP, nothing, debug, verb) #! ~45GB
    GC.gc(); solve_alm = @elapsed res, _, _ = alm(x0, lambda0, mux0, sigma0, tau0, PI, AP, SP, nothing, debug, verb, fullit) #! ~45GB
    writeout_alm(res, RESULTPATH, PI, outname, [m], n, 1)
  end
  # solve grb
  begin
    GC.gc(); solve_grb = @elapsed res = solve_cvar_reduced_gurobi(PI, options) #! >3600s (5285s), ~80GB
    res[:y] = vcat([PI.A[l] * res[:x] .+ PI.b[l] for l in 1:PI.L]...)
    calcs_manual!(
      res, obj, PI,
      res[:x], res[:y], res[:z], res[:t],
      res[:lambda], res[:mu]
    )
    writeout_grb(res, RESULTPATH * "/cold/", PI, outname, [m], n, 1)
  end
  # solve grb_qr
  begin
    # GC.gc(); solve_grb_qr = @elapsed res = solve_qr_gurobi(-A, yy, k, options) #! >3600s, ~120GB
    GC.gc(); solve_grb_qr = @elapsed res = solve_qr_gurobi(PI, options) #! >3600s, ~120GB
    writeout_grb_qr(res, RESULTPATH * "/cold/", PI, outname, [m], n, 1)
  end
  # solve grb_oa; needs to be last bc it modifies PI.xlo & PI.xhi
  begin
    PI.xlo .= -10.0
    PI.xhi .= +10.0
    GC.gc(); solve_grb_oa = @elapsed res = solve_cvar_oa_reduced_gurobi(PI, options) #! >3600s, ~80GB
    res[:y] = vcat([PI.A[l] * res[:x] .+ PI.b[l] for l in 1:PI.L]...)
    calcs_manual!(
      res, obj, PI,
      res[:x], res[:y], res[:z], res[:t],
      res[:lambda], res[:mu]
    )
    writeout_grb_oa(res, RESULTPATH * "/cold/", PI, outname, [m], n, 1)
  end
end

#
# path timing
#

global const TOL = 1e-4
r = 0.0
k = 0
PI = ProblemInstance(
  -1, -1, n, [m], 1, [A], [b],
  [k], #! overwrite
  [r],
  Function[Atmatinv_action(m, k)], #! overwrite
  fill(-Inf, n), fill(+Inf, n), fill(false, n), fill(false, n), true,
  1, nothing, nothing, c[1:n], Dict()
)

Tf = Float64
obj = FunctionValue(Tf, PI.n, PI.objtype) # 2==diagonal
obj.g .= PI.c
BLAS.set_num_threads(nthreads)
x0 = zeros(Tf, PI.n);
x0[end] = 1.0;
lambda0 = [zeros(m[l]) for l in eachindex(m)];
mux0 = zeros(Tf, PI.n);
debug = false
verb = 2
sigma0 = 1e-5 #! worse: sigma0 = 1e-8
tau0 = 1e-6   #! worse: tau0 = 1e-9

AP.stoptol = 1e-4
SP.maxlit = 20
SP.ls_method=:sw
pctgrid = [1e-3, 1e-2, collect(range(0.025, stop=0.975, step=0.025))..., 1-1e-2, 1-1e-3]
kgrid = [[Int(ceil(pct*PI.m[1]))] for pct in pctgrid]
solve_alm_grid = @elapsed res_alm_grid_full = alm(x0, lambda0, mux0, sigma0, tau0, PI, AP, SP, kgrid, debug, verb) #! ~45GB
writeout_alm_grid(res_alm_grid_full, RESULTPATH * "/warm/", PI, [m], n, 1)
zero_x = true
solve_alm_grid = @elapsed res_alm_grid_full_zero = alm(x0, lambda0, mux0, sigma0, tau0, PI, AP, SP, kgrid, debug, verb, false, zero_x) #! ~45GB

solve_alm_grid_final = @elapsed res_alm_grid_final = alm(x0, lambda0, mux0, sigma0, tau0, PI, AP, SP, kgrid[41:43], debug, verb) #! ~45GB
solve_alm_grid_initial = @elapsed res_alm_grid_initial = alm(x0, lambda0, mux0, sigma0, tau0, PI, AP, SP, kgrid[1:2], debug, verb) #! ~45GB

solve_alm_oa = @elapsed res_alm_oa = alm_oa(x0, lambda0, mux0, sigma0, tau0, PI, AP, SP, 1.0, 0.5, debug, verb) #! ~45GB

ALM   GRB   GRB-OA PSG  OSQP
3.7e1 3.3e2 2.6e1 7.7e2 3.6e3 
2.0e0 1.6e1 2.3e0 3.4e1 3.6e2 
3.6e1 1.2e2 8.4e1 3.6e3 3.6e3

1.2e2 1.8e2 2.1e2 6.8e2 3.6e3 
4.1e0 1.3e1 1.5e1 6.2e1 3.5e2 
6.5e1 1.2e2 2.8e2 5.7e2 2.5e3

5.9e0 1.5e2 1.1e1 4.7e0 5.8e2 
2.8e-1 5.8e0 1.2e0 2.2e-1 1.4e2 
2.1e0 2.6e1 2.6e0 3.3e1 4.1e2

1.2e1 1.3e2 1.9e2 8.1e0 8.8e2 
3.0e-1 4.3e0 5.4e0 4.4e-1 1.0e2 
2.8e0 2.1e1 4.7e1 5.6e1 5.4e2

ALM   GRB   GRB-OA PSG  OSQP
3.5e1 5.3e2 3.2e1 9.9e0 3.6e3 
1.5e0 2.0e1 2.8e0 6.7e-1 1.9e2 
6.6e1 1.3e2 1.2e2 2.3e1 1.8e3

6.4e1 2.0e2 3.3e2 1.3e1 3.6e3 
3.6e0 1.5e1 2.0e1 8.2e-1 4.0e2 
1.0e2 1.3e2 2.7e2 1.7e1 2.3e3

2.1e0 2.3e2 1.5e1 2.1e0 2.4e2 
1.6e-1 8.7e0 9.9e-1 1.2e-1 3.7e1 
1.5e0 2.5e1 3.2e0 7.5e-1 2.1e2

3.1e0 1.4e2 1.6e2 2.5e0 1.1e3 
1.6e-1 5.8e0 6.0e0 1.3e-1 2.6e1 
2.8e0 2.4e1 4.7e1 1.0e0 2.8e2
