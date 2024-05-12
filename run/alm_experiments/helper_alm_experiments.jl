#
# helpers
#

mutable struct TestInstance{Tf<:AbstractFloat,Vf<:AbstractVector{Tf},Ti<:Integer,Vi<:AbstractVector{Ti}}
  seed::Ti
  objtype::Ti
  n::Ti
  m::Vi
  pct::Vf
  scale_col::Bool
  scale_row::Bool
  t::Tf # mks constraint level: 0 = hard, 1 = easy
  
  sigma0::Tf
  tau0::Tf
  AP::ALMParameters
  SP::SSNParameters
end
get_k(TI::TestInstance) = Int.(ceil.(TI.pct .* TI.m))

function antidiagonal(base::Integer, lo::Integer, hi::Integer, step::Integer=1)
  #=
   1  1  1
   1  1  0
   1  0  0
  =#
  lo_hi = lo:step:hi
  n = length(lo_hi)
  idx = zeros(Bool, n, n)
  val = Array{Array}(undef, n, n)
  for i in 1:n
    for j in 1:n
      if i+j <= length(lo_hi)+1
        val[i,j] = [base^lo_hi[i], base^lo_hi[j]]
        idx[i,j] = true
      end
    end
  end
  return val[idx], idx
end
function antidiagonal(
  base_m::Integer, lo_m::Integer, hi_m::Integer, step_m::Integer,
  base_n::Integer, lo_n::Integer, hi_n::Integer, step_n::Integer  
)
  #=
   1  1  1
   1  1  0
   1  0  0
  =#
  lo_hi_m = lo_m:step_m:hi_m
  lo_hi_n = lo_n:step_n:hi_n
  nrow = length(lo_hi_m)
  ncol = length(lo_hi_n)
  @assert(nrow == ncol)
  ndim = nrow
  idx = zeros(Bool, nrow, ncol)
  val = Array{Array}(undef, nrow, ncol)
  for i in 1:nrow
    for j in 1:ncol
      if i+j <= ndim+1
        val[i,j] = [base_m^lo_hi_m[i], base_n^lo_hi_n[j]]
        idx[i,j] = true
      end
    end
  end
  return val[idx], val, idx
end

function get_dim(nterms::Integer)
  return Int(0.5 * (sqrt(8*nterms+1)-1)) # solve n(n+1)/2 = length(mn)
end

function DelimitedFiles.writedlm(TI::TestInstance, fname::String)
  data = OrderedDict(
    [k => getfield(TI, k) for k in fieldnames(TestInstance)[1:end-2]]...,
    [k => getfield(TI.AP, k) for k in fieldnames(ALMParameters)]...,
    [k => getfield(TI.SP, k) for k in fieldnames(SSNParameters)]...,
  )
  data[:lin_solver]=data[:lin_solver][1]
  writedlm(fname, data, '|')
end
function reconstruct_TestInstance(fname::String)
  idx = 0
  x = readdlm(fname, '|')
  x[16,2] = Symbol(x[16,2])
  x[21,2] = Symbol(x[21,2])
  idx += 1
  seed = x[idx,2]
  idx += 1
  objtype = x[idx,2]
  idx += 1
  n = x[idx,2]
  idx += 1
  m = parse.(Int64,split(x[idx,2], ('[',']',','))[2:end-1])
  idx += 1
  pct = parse.(Float64,split(x[idx,2], ('[',']',','))[2:end-1])
  idx += 1
  scale_col = x[idx,2]
  idx += 1
  scale_row = x[idx,2]
  idx += 1
  t = x[idx,2]
  idx += 1
  sigma0 = x[idx,2]
  idx += 1
  tau0 = x[idx,2]
  idx += 1
  AP = ALMParameters(x[idx,2], x[idx+1,2])
  idx += 2
  SP = SSNParameters(x[idx:idx+7,2]..., (x[idx+8,2],1))
  return TestInstance(
    seed, objtype, n, m, pct, scale_col, scale_row, t, sigma0, tau0,
    AP, SP
  )
end

function writeout_alm(res::Dict, nactive::Integer, datapath::String, PI::ProblemInstance, outname::String, m::Vector, n::Integer, repi::Integer)
  # common
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__pinfeas.csv", res[:pinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__dinfeas.csv", res[:dinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__relgap.csv", res[:relgap])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__walltime.csv", res[:walltime])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__walltime_i1.csv", haskey(res, :walltime_i1) ? res[:walltime_i1] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__walltime_i2.csv", haskey(res, :walltime_i2) ? res[:walltime_i2] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__pobj.csv", res[:pobj])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__dobj.csv", res[:dobj])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__iter.csv", res[:iter])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__status.csv", res[:status])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__avg_step_size.csv", res[:avg_step_size])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__nnz.csv", sum(prod.(size.(PI.A))))
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__nactive.csv", nactive)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__kktres.csv", res[:err])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__r.csv", res[:maxksumrhs])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__mks.csv", res[:maxksum])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__cvar.csv", res[:cvar]) # vector
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__cvar_infeas.csv", res[:cvar_infeas]) # scalar
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__cvar_rhs.csv", res[:cvar_rhs]) # vector
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__normx.csv", res[:normx]) # scalar
  # specific
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__retcode.csv", res[:retcode])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__ssnstatus.csv", res[:substatus])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__feval.csv", sum(res[:lsit]))
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__vphieval.csv", sum(res[:lsit]))
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__gphieval.csv", sum(res[:nit]))
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__Hphieval.csv", sum(res[:nit]))
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__avg_nd_size.csv", res[:avg_nd_size])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__absbeta.csv", res[:absbeta])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__absalpha.csv", res[:absalpha])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__ssnpinfeas.csv", res[:sub_pinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__ssndinfeas.csv", res[:sub_dinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__ssnrelgap.csv", res[:sub_relgap])
  # timing
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_proj.csv", res[:ctime_proj])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_sort.csv", res[:ctime_sort])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_nd.csv", res[:ctime_nd])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_Hphi.csv", res[:ctime_Hphi])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_gphi.csv", res[:ctime_gphi])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_vphi.csv", res[:ctime_vphi])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_dirn.csv", res[:ctime_nd] + res[:ctime_Hphi]) #! combined
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_A.csv", res[:ctime_Ax])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_At.csv", res[:ctime_Atlam])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_ssn.csv", res[:timessn])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__lsit.csv", res[:lsit])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__sps.csv", res[:sps])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__psc.csv", res[:partial_sort_count])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__nit.csv", res[:nit])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__sigma.csv", res[:sigma])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__tau.csv", res[:tau])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__pct_timed.csv", res[:pct_timed])
end
function writeout_grb(res::Dict, nactive::Integer, datapath::String, PI::ProblemInstance, outname::String, m::Vector, n::Integer, repi::Integer, method::String="")
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__pinfeas.csv", res[:pinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__dinfeas.csv", res[:dinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__relgap.csv", res[:relgap])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__walltime.csv", res[:walltime])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__pobj.csv", res[:pobj])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__iter.csv", res[:iter])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__status.csv", res[:status])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__t.csv", res[:t])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__nnz.csv", res[:nnz])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__nactive.csv", nactive)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__mks.csv", res[:maxksum])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__r.csv", res[:maxksumrhs])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__cvar.csv", res[:cvar]) # vector
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__cvar_infeas.csv", res[:cvar_infeas]) # scalar
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__cvar_rhs.csv", res[:cvar_rhs]) # vector
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb$(method)__normx.csv", res[:normx]) # scalar
end
function writeout_grb_oa(res::Dict, nactive::Integer, datapath::String, PI::ProblemInstance, outname::String, m::Vector, n::Integer, repi::Integer, method::String="")
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__pinfeas.csv", res[:pinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__dinfeas.csv", res[:dinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__relgap.csv", res[:relgap])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__walltime.csv", res[:walltime])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__subtimes.csv", res[:times])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__subiters.csv", res[:iters])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__pobj.csv", res[:pobj])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__iter.csv", res[:iter])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__status.csv", res[:status])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__ncon.csv", res[:numcon])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__mks.csv", res[:maxksum])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__r.csv", res[:maxksumrhs])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__cvar.csv", res[:cvar]) # vector
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__cvar_infeas.csv", res[:cvar_infeas]) # scalar
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__cvar_rhs.csv", res[:cvar_rhs]) # vector
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa$(method)__normx.csv", res[:normx]) # scalar
end
# function writeout_grb_benchmark(res::Dict, nactive::Integer, datapath::String, PI::ProblemInstance, outname::String, m::Vector, n::Integer, repi::Integer)
#   writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_benchmark__walltime.csv", res[:walltime])
#   writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_benchmark__pobj.csv", res[:pobj])
#   writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_benchmark__iter.csv", res[:iter])
#   writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_benchmark__status.csv", res[:status])
#   writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_benchmark__t.csv", res[:t])
#   writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_benchmark__nactive.csv", nactive)
# end
function writeout_osqp(res::Dict, nactive::Integer, datapath::String, PI::ProblemInstance, outname::String, m::Vector, n::Integer, repi::Integer)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__pinfeas.csv", res[:pinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__dinfeas.csv", res[:dinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__relgap.csv", res[:relgap])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__walltime.csv", res[:walltime])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__pobj.csv", res[:pobj])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__iter.csv", res[:iter])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__status.csv", res[:status])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__nnz.csv", res[:nnz])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__nactive.csv", nactive)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__mks.csv", res[:maxksum])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__r.csv", res[:maxksumrhs])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__cvar.csv", res[:cvar]) # vector
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__cvar_infeas.csv", res[:cvar_infeas]) # scalar
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__cvar_rhs.csv", res[:cvar_rhs]) # vector
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__osqp__normx.csv", res[:normx]) # scalar
end

function generate_problem(TI::TestInstance, repi, experi, Tf=Float64)
  # PI = nothing; CI = nothing; OI = nothing; GC.gc()
  prob_time = @elapsed begin
  seed = TI.seed
  CI = ConstraintInstance(seed, TI.n, TI.m, get_k(TI), TI.scale_row, TI.scale_col, TI.t) #! get rid of TI.t?
  OI = ObjectiveInstance(seed, TI.objtype, TI.n)
  PI = ProblemInstance(CI, OI)
  println("============================================================================================================")
  println("objtype = $(PI.objtype); seed = $seed; experi = $experi; L = $(length(TI.m)); pct = $(TI.pct); k = $(get_k(TI)); m = $(TI.m[1]); n = $(TI.n[1])")
  println("============================================================================================================")

  # obj function values (phi is initialized inside alm)
  obj = FunctionValue(Tf, TI.n, TI.objtype) # 2==diagonal
  if !isnothing(PI.C)
    obj.H .= PI.C
  else
    obj.g .= PI.c
  end; end # prob time
  println("problem generation time: $(round(prob_time, digits=2))s")
  return PI, obj
end

function test_alm_timing(
  pct::Float64, L::Int64, mtd::Vector{Symbol}, datapath::String,
  mns::Vector, nrep::Integer,
  TI::TestInstance, solver_options::Dict, Tf::DataType, writeoutput::Bool=true,
)
  # initialize
  @assert(pct <= 1.0 && pct >= 0.0)
  outname = "k_$(Int(pct*100))pct__L_$L"
  if !isdir(datapath * "obj_$(TI.objtype)/$outname/"); mkpath(datapath * "obj_$(TI.objtype)/$outname/"); end
  # psg_solvers = ["VAN", "CAR", "BULDOZER", "TANK", "VANGRB", "CARGRB", "HELI"]
  psg_solvers = ["VANGRB"] #! GRB ones cannot handle >= 1000 variables? VANGRB fails for m=16k, n=500
  psg_solvers = ["HELI"] #! VANGRB fails for m=16k, n=500
  osqp_cutoff = 10_000_000
  mks_active_tol = 1e-5
  maxtime_original = copy(options[:maxtime])
  maxnit0 = TI.SP.maxnit
  stagit0 = TI.SP.stagit
    
  # rep function
  sim_repi(repi::Int64) = begin
    
    # solve problems
    for experi in eachindex(mns)

      #
      # setup experiment
      #

      # experiment data
      m = mns[experi][1]
      n = mns[experi][2]
      TI.m = [m for l in 1:L]
      TI.n = n

      # reset cvar-solver options
      options = copy(solver_options);

      #! skip cases
      # if TI.objtype==2
      #   if !(pct==0.01 && L==1)
      #     continue
      #   end
      # end
      # if !(pct==0.01 && L==1 && TI.objtype==1 && TI.n==1000 && TI.m[1]==100_000)
      #   continue
      # end
      # if pct != 0.01 && L != 16
      #   continue
      # end
      # if TI.objtype == 2
      #   continue
      # end
      # if experi <= 12
      #   continue
      # end

      #
      # solve alm (hp = highprecision)
      #
      
      if :alm ∈ mtd
        # generate problem
        PI = nothing; GC.gc()
        PI, obj = generate_problem(TI, repi, experi, Float64)

        BLAS.set_num_threads(nthreads) # global const!
        x0 = zeros(Tf, TI.n);
        # x0 = -PI.Cinv * PI.c
        # x0 = rand(Tf, TI.n)
        lambda0 = [zeros(TI.m[l]) for l in eachindex(TI.m)];
        mux0 = zeros(Tf, TI.n);
        debug = false
        verb = 1
        TI.tau0 = 1.0
        kgrid = nothing
        fullit = sum(PI.m)/sum(PI.L) <= 2PI.n
        TI.SP.maxnit = maxnit0
        TI.SP.stagit = sum(PI.m)/PI.L > 4PI.n ? stagit0 : 100
        AP.BXratio = sqrt(sum(PI.m)/sum(PI.L)/PI.n)
        sigma0 = 1 / AP.BXratio
        tau0 = sigma0 / AP.BXratio
        solve_alm = @elapsed res, _, _ = alm(x0, lambda0, mux0, sigma0, tau0, PI, TI.AP, TI.SP, nothing, debug, verb, fullit)
        nactive = sum(res[:maxksumrhs] - res[:maxksum] .<= mks_active_tol)
        println("solve time alm: $(round(solve_alm, digits=2))s")
        println("alm solved? $(res[:retcode])")

        # write output
        if writeoutput
          writeout_alm(res, nactive, datapath, PI, outname, TI.m, TI.n, repi)
        end
        partition = get_partition(res[:y], PI.k)
        PI = nothing; GC.gc()
      end

      #
      # solve gurobi
      #

      # if !(m >= 100_000 && n >= 10_000)
      if :grb ∈ mtd
        # generate problem
        PI = nothing; GC.gc()
        PI, obj = generate_problem(TI, repi, experi, Float64)
        if PI.objtype == 1 || PI.objtype == 2
          # options[:method] = -1 # slow/memory issues
          # options[:method] = 2 # 0 = primal simplex, 1 = dual simplex, 2 = barrier
          options[:method] = 1 # 0 = primal simplex, 1 = dual simplex, 2 = barrier
          options[:crossover] = 0
          options[:nthreads] = 1 # need two: one for simplex and one for IPM; can't handle large problems though
        end
        solve_grb = @elapsed res = solve_cvar_reduced_gurobi(PI, options, true)
        println("solve time gurobi: $(round(solve_grb, digits=2))s")
        if haskey(res, :x)
          PI, obj = generate_problem(TI, repi, experi, Float64)
          res[:y] = vcat([PI.A[l] * res[:x] .+ PI.b[l] for l in 1:L]...)
          calcs_manual!(
            res, obj, PI,
            res[:x], res[:y], res[:z], res[:t],
            res[:lambda], res[:mu]
          )
          nactive = sum(res[:maxksumrhs] - res[:maxksum] .<= mks_active_tol)
          
          # write output
          if writeoutput
            writeout_grb(res, nactive, datapath, PI, outname, TI.m, TI.n, repi, (options[:method] == 1 ? "ds" : ""))
          end
          PI = nothing; GC.gc()
        end
      end

        if :grb_oa ∈ mtd
        # generate problem
        PI = nothing; GC.gc()
        PI, obj = generate_problem(TI, repi, experi, Float64)
        if PI.objtype == 1 || PI.objtype == 2
          # options[:method] = -1 # slow/memory issues
          # options[:method] = 2 # 0 = primal simplex, 1 = dual simplex, 2 = barrier
          options[:method] = 1 # 0 = primal simplex, 1 = dual simplex, 2 = barrier
          options[:crossover] = 0
          options[:nthreads] = 1 # need two: one for simplex and one for IPM; can't handle large problems though
        end
        solve_grb_oa = @elapsed res = solve_cvar_oa_reduced_gurobi(PI, options)
        println("solve time gurobi-oa: $(round(solve_grb_oa, digits=2))s")
        if haskey(res, :x)
          PI, obj = generate_problem(TI, repi, experi, Float64)
          res[:y] = vcat([PI.A[l] * res[:x] .+ PI.b[l] for l in 1:L]...)
          calcs_manual!(
            res, obj, PI,
            res[:x], res[:y], res[:z], res[:t],
            res[:lambda], res[:mu]
          )
          nactive = sum(res[:maxksumrhs] - res[:maxksum] .<= mks_active_tol)
          
          # write output
          if writeoutput
            writeout_grb_oa(res, nactive, datapath, PI, outname, TI.m, TI.n, repi, (options[:method] == 1 ? "ds" : ""))
          end
          PI = nothing; GC.gc()
        end
      end
      # end

      #
      # solve PSG
      #

      if Sys.iswindows() #! PSG only works on windows
        GC.gc()
        if :psg ∈ mtd
          for psg_solver in psg_solvers
            PI = nothing; GC.gc()
            PI, obj = generate_problem(TI, repi, experi, Float64)

            # clean up temporary directories
            for fi in readdir(psg_problemdata)
              rm(psg_problemdata * fi)
            end
            for fi in readdir(psg_problemoutput)
              rm(psg_problemoutput * fi)
            end

            # transfer constraint data
            writedlm(psg_problemdata * "L.csv", PI.L)
            for l in 1:L
              writedlm(psg_problemdata * "A$(l).csv", PI.A[l])
              writedlm(psg_problemdata * "b$(l).csv", PI.b[l])
            end
            writedlm(psg_problemdata * "k.csv", PI.k)
            writedlm(psg_problemdata * "r.csv", PI.r)
            writedlm(psg_problemdata * "lb.csv", PI.xlo)
            writedlm(psg_problemdata * "ub.csv", PI.xhi)

            # transfer objective data
            if !isnothing(PI.C)
              writedlm(psg_problemdata * "Cq_diag.csv", diag(PI.C))
              writedlm(psg_problemdata * "Cq.csv", 0.0)
            else
              writedlm(psg_problemdata * "Cq_diag.csv", zeros(PI.n))
            end
            writedlm(psg_problemdata * "cl.csv", PI.c)

            # transfer solver
            if PI.n > 10_000
              writedlm(psg_problemdata * "SOLVER.csv", ["CARGRB"])
            else
              writedlm(psg_problemdata * "SOLVER.csv", [psg_solver])
            end
            # if PI.objtype == 1
            #   writedlm(psg_problemdata * "SOLVER.csv", ["CARGRB"]) # http://www.aorda.com/html/PSG_Help_HTML/index.html?solver.htm
            # end

            # call matlab to solve psg
            # 2021b: C:\Program Files\MATLAB\R2021b\bin
            # from command prompt: "C:\Program Files\MATLAB\R2021b\bin\matlab.exe" -nosplash -nodesktop  -r "run C:\Users\roth0674\Documents\GitHub\cc-non\code\run\alm_experiments\run_psg.m" 
            # run(`matlab -nosplash -nodesktop -r "run 'C:\Users\roth0674\Documents\GitHub\cc-non\code\run\alm_experiments\test_psg.m'"`)
            # close(MATLAB_session)
            psg_observedtime = @elapsed begin
              # MATLAB_session = MSession()
              eval_string(MATLAB_session,"run run_psg.m")
            end

            # get stats
            psg_xbar = vec(readdlm(psg_problemoutput * "xbar.csv",','))
            psg_y = [PI.A[l] * psg_xbar .+ PI.b[l] for l in 1:PI.L]
            psg_cvar_rhs = vec(PI.r ./ PI.k)
            
            # internal
            psg_cvar_internal = vec(readdlm(psg_problemoutput * "mks.csv",',')) # is value of constraint (ie cvar) not mks
            psg_mks_internal = psg_cvar_internal .* PI.k
            psg_cvar_infeas_internal = maximum(max.(psg_cvar_internal .- psg_cvar_rhs, 0.0))
            psg_mks_infeas_internal = maximum(max.(psg_mks_internal .- PI.r, 0.0))
            psg_pinfeas_internal = max(
              psg_cvar_infeas_internal, # cvar
              maximum(max.(PI.xlo .- psg_xbar, 0.0)), # box
              maximum(max.(psg_xbar .- PI.xhi, 0.0)), # box
            )
            
            # computed
            psg_mks = [maxksum(psg_y[l], PI.k[l]) for l in 1:PI.L]
            psg_cvar = psg_mks ./ PI.k
            psg_cvar_infeas = maximum(max.(psg_cvar .- psg_cvar_rhs, 0.0))
            psg_mks_infeas = maximum(max.(psg_mks .- PI.r, 0.0))
            psg_pinfeas = max(
              psg_cvar_infeas, # cvar
              maximum(max.(PI.xlo .- psg_xbar, 0.0)), # box
              maximum(max.(psg_xbar .- PI.xhi, 0.0)), # box
            )
            if PI.objtype == 2
              psg_pobj = 0.5 * psg_xbar' * PI.C * psg_xbar + PI.c' * psg_xbar
            else
              psg_pobj = PI.c' * psg_xbar
            end
            psg_pobj_internal = vec(readdlm(psg_problemoutput * "pobj.csv",','))[1]
            psg_gap = vec(readdlm(psg_problemoutput * "gap.csv",','))[1]
            psg_load_time = vec(readdlm(psg_problemoutput * "load_time.csv",','))[1]
            psg_preprocess_time = vec(readdlm(psg_problemoutput * "preprocess_time.csv",','))[1]
            psg_solve_time = vec(readdlm(psg_problemoutput * "solve_time.csv",','))[1]
            psg_status = vec(readdlm(psg_problemoutput * "status.csv",','))[1]
            psg_nactive = sum(PI.r - psg_mks .<= mks_active_tol)

            # writedata
            if writeoutput
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__observedtime.csv", psg_observedtime)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__mks_internal.csv", psg_mks_internal)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__cvar_internal.csv", psg_cvar_internal)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__mks_infeas_internal.csv", psg_mks_infeas_internal)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__cvar_infeas_internal.csv", psg_cvar_infeas_internal)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__pinfeas_internal.csv", psg_pinfeas_internal)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__mks.csv", psg_mks)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__cvar.csv", psg_cvar)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__mks_infeas.csv", psg_cvar_infeas)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__cvar_infeas.csv", psg_cvar_infeas)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__pinfeas.csv", psg_pinfeas)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__r.csv", PI.r)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__cvar_rhs.csv", psg_cvar_rhs)
              # writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__xbar.csv", psg_xbar)
              # writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__dinfeas.csv", res[:dinfeas]) #! NA
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__relgap.csv", psg_gap)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__walltime.csv", psg_preprocess_time + psg_solve_time)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__solvetime.csv", psg_solve_time)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__pobj.csv", psg_pobj)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__pobj_internal.csv", psg_pobj_internal)
              # writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg__iter.csv", res[:iter]) #! NA
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__status.csv", psg_status == "optimal" && psg_pinfeas_internal <= options[:pinfeas_tol_hi])
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__status_internal.csv", [psg_status])
              # writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__nnz.csv", -1) #! NA
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__nactive.csv", psg_nactive)
              writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__psg$(psg_solver)__normx.csv", norm(psg_xbar))
            end
            PI = nothing; GC.gc()
          end
        end
      end

      #
      # solve osqp
      #
      
      # if :osqp ∈ mtd && TI.n * sum(TI.m) <= osqp_cutoff
      if :osqp ∈ mtd
        # generate problem
        PI, obj = generate_problem(TI, repi, experi, Float64)

        options[:lin_sys_solver] = "qdldl"
        options[:maxtime] = copy(maxtime_original)
        solve_osqp = @elapsed res = solve_cvar_reduced_osqp(PI, options)
        # solve_osqp = @elapsed res_reduced = solve_cvar_reduced_osqp(PI, options)
        # calcs_manual!(
        #   res_reduced, obj, PI,
        #   res_reduced[:x], res_reduced[:y], res_reduced[:z], res_reduced[:t],
        #   res_reduced[:lambda], res_reduced[:mu]
        # )
        println("solve time osqp: $(round(solve_osqp, digits=2))s")
        if haskey(res, :x)
          PI, obj = generate_problem(TI, repi, experi, Float64)
          res[:y] = vcat([PI.A[l] * res[:x] .+ PI.b[l] for l in 1:L]...)
          calcs_manual!(
            res, obj, PI,
            res[:x], vcat(res[:y]...), vcat(res[:z]...), res[:t],
            vcat(res[:lambda]...), res[:mu]
          )
          nactive = sum(res[:maxksumrhs] - res[:maxksum] .<= mks_active_tol)
          
          # write output
          writeout_osqp(res, nactive, datapath, PI, outname, TI.m, TI.n, repi)
          PI = nothing; GC.gc()
        end
      end
    end # exper_i (n,m)
  end # sim_repi
  
  # run experiments
  if Sys.iswindows()
    pmap(sim_repi, 1:nrep)
  else
    sim_repi(nrep)
  end
  
  # return
  return outname
end

function load_alm_results(datapath::String, outname::String)
  out = Dict()
  expers = readdir(datapath * "/" * outname * "/")
  expers = expers[[e[1] .!= '.' for e in expers]]
  mnr = [x[1] for x in split.(expers,"__")]
  repid = [parse(Int64, x[3]) for x in split.(mnr,"_")]
  nrep = maximum([parse(Int64, x[3]) for x in split.(mnr,"_")])
  mtd = [x[2] for x in split.(expers,"__")]
  mn = [split(x[1], "_")[1] * "_" * split(x[1], "_")[2] for x in split.(expers,"__")]
  fld = [replace(x[3], ".csv"=>"") for x in split.(expers,"__")]
  
  mtd_names = unique(mtd)
  mnr_names = unique(mnr)
  mn_names = unique([x[1]*"_"*x[2] for x in split.(mnr, "_")])

  # initialize dictionaries
  for i in eachindex(mtd_names)
    out[mtd_names[i]] = Dict()
    for j in eachindex(mn_names)
      out[mtd_names[i]][mn_names[j]] = Dict()
    end
  end

  # populate dictionaries
  for i in eachindex(expers)
    mtd_name = mtd[i]
    mn_name = mn[i]
    rep_name = repid[i]
    fld_name = fld[i]
    file_name = expers[i]
    rep_file_name_mask = expers .∈ Ref([mn_name * "_" * string(repi) * "__" * string(mtd_name) * "__" * fld_name * ".csv" for repi in 1:nrep])
    data = []
    println(expers[i])
    for (repi,fn) in enumerate(expers[rep_file_name_mask])
      println(fn)
      rec = readdlm(datapath * "/" * outname * "/" * fn)
      if eltype(rec[1])==Char
        if repi == 1; data = String[]; end
        push!(data, rec[1])
      elseif length(rec) == 1
        if repi == 1; data = Float64[]; end
        push!(data, rec[1])
      else
        if repi == 1; data = Vector{Float64}[]; end
        push!(data, vec(rec))
      end
    end
    out[mtd_name][mn_name][fld_name] = data
  end
  return out
end