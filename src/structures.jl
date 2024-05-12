#=
  data for primal problem instance:
  min_x { f(x) : y = Ax + b, M_k(y) <= r }
=#
mutable struct ConstraintInstance{
  Tf<:AbstractFloat, Ti<:Integer,
  Mf<:AbstractMatrix{Tf}, Vf<:AbstractVector{Tf},
  VMfF<:Union{Vector{Mf},Vector{Function}}
}
  # problem id
  seed  ::  Ti

  # problem dimensions
  n  ::  Ti # x-dim
  m  ::  Vector{Ti} # y-dim
  L  ::  Ti # number of CVaR constraints
  
  # CVaR constraints
  A  ::  Vector{Mf} # m × n
  b  ::  Vector{Vf} # m
  k  ::  Vector{Ti} # threshold
  r  ::  Vector{Tf} # rhs value
  
  # box constraints
  Bt  ::  VMfF # m × m = Atildeinv
  xlo  ::  Vf # xlo <= x
  xhi  ::  Vf # x <= xhi
  xloidx  :: Vector{Bool}
  xhiidx  :: Vector{Bool}
  Xunconstr  ::  Bool
end
mutable struct ObjectiveInstance{
  Tf<:AbstractFloat, Ti<:Integer,
  Vf<:AbstractVector{Tf}, NDf<:Union{Nothing,Diagonal{Tf}}
}
  seed  ::  Ti
  objtype  ::  Ti # objective type: 1 = linear, 2 = quadratic (pd diagonal), 3 = quadratic (pd), 4 = general
  C  ::  NDf
  Cinv  ::  NDf
  c  ::  Vf
end
struct ProblemInstance{
  Tf<:AbstractFloat, Ti<:Integer,
  Mf<:AbstractMatrix{Tf}, Vf<:AbstractVector{Tf}, NDf<:Union{Nothing,Diagonal{Tf}},
  VMfF<:Union{Vector{Mf},Vector{Function}}
}
  # problem id
  conseed  ::  Ti
  objseed  ::  Ti
  
  # problem dimensions
  n  ::  Ti # x-dim
  m  ::  Vector{Ti} # y-dim
  L  ::  Ti # number of CVaR constraints
  
  # CVaR constraints
  A  ::  Vector{Mf} # m × n
  b  ::  Vector{Vf} # m
  k  ::  Vector{Ti} # threshold
  r  ::  Vector{Tf} # rhs value
  
  # box constraints
  Bt  ::  VMfF # m × m = Atildeinv
  xlo  ::  Vf # xlo <= x
  xhi  ::  Vf # x <= xhi
  xloidx  :: Vector{Bool}
  xhiidx  :: Vector{Bool}
  Xunconstr  ::  Bool

  # objective
  objtype  ::  Ti # objective type: 1 = linear, 2 = quadratic (pd diagonal), 3 = quadratic (pd), 4 = general
  C  ::  NDf
  Cinv  ::  NDf
  c  ::  Vf

  # scaling
  scale::Dict
end
function ProblemInstance(CI::ConstraintInstance, OI::ObjectiveInstance, Tfr::DataType=Float64)
  L = length(CI.A)
  return ProblemInstance(
    CI.seed, OI.seed,
    CI.n, CI.m, CI.L, [Tfr.(CI.A[l]) for l in 1:L], [Tfr.(CI.b[l]) for l in 1:L], CI.k, Tfr.(CI.r), CI.Bt, Tfr.(CI.xlo), Tfr.(CI.xhi), CI.xloidx, CI.xhiidx, CI.Xunconstr,
    OI.objtype, (!isnothing(OI.C) ? Tfr.(OI.C) : OI.C), (!isnothing(OI.Cinv) ? Tfr.(OI.Cinv) : OI.Cinv), Tfr.(OI.c), Dict()
  )
end
function rescale!(PI::ProblemInstance, scale_k::Bool, scale_ruiz::Bool, ruiz_maxit::Integer=25, ruiz_tol::Real=1e-6)
  if scale_ruiz
    A = vcat(PI.A...)
    invE_, _, invD = ruiz!(A, ruiz_maxit, ruiz_tol)
    invE = [zeros(PI.m[l]) for l in 1:PI.L]
    for l in 1:PI.L
      idx = ([1;1 .+ cumsum(PI.m[1:end-1])], cumsum(PI.m))
      PI.A[l] .= A[idx[1][l]:idx[2][l],:]
      invE[l] .= invE_[idx[1][l]:idx[2][l]]
    end
    #! separate scaling not correct
    # scaling = [ruiz!(PI.A[l], ruiz_maxit, ruiz_tol) for l in 1:PI.L]
    # invE = [scaling[l][1] for l in 1:PI.L]
    # invD = [scaling[l][3] for l in 1:PI.L]
  end
  if scale_k
    for l in 1:PI.L
      PI.A[l] .*= 1/PI.k[l]
      PI.b[l] .*= 1/PI.k[l]
      PI.r[l] *= 1/PI.k[l]
    end
  end
  if scale_ruiz
    PI.scale[:invE] = invE
    PI.scale[:invD] = invD
  end
end

#=
  container to hold function's 0th (v), 1st (g), and 2nd (H) order information
=#
mutable struct FunctionValue{Tf<:AbstractFloat, Vf<:AbstractVector{Tf}, Mf<:AbstractMatrix{Tf}}
  v :: Tf
  g :: Vf
  H :: Mf
  d :: Tf
end
function FunctionValue(Tf::DataType, n::Integer, objtype::Integer)
  if objtype == 1
    # return FunctionValue(zero(Tf), zeros(Tf, n), spzeros(Tf, n, n), zero(Tf))
    return FunctionValue(zero(Tf), zeros(Tf, n), Diagonal(zeros(Tf, n)), zero(Tf))
  elseif objtype == 2
    return FunctionValue(zero(Tf), zeros(Tf, n), Diagonal(zeros(Tf, n)), zero(Tf))
  else
    return FunctionValue(zero(Tf), zeros(Tf, n), zeros(Tf, n, n), zero(Tf))
  end
end

#=
  SSN and ALM parameters
=#
mutable struct SSNParameters{Tf<:AbstractFloat, Ti<:Integer}
  maxnit  ::  Ti
  maxlit  ::  Ti
  tolsub  ::  Tf
  ls_method  ::  Symbol
  ls_tau  ::  Tf
  ls_c1  ::  Tf
  ls_c2  ::  Tf
  stagit  ::  Ti # number of iteration until stagnate
  lin_solver  ::  Tuple
end
function SSNParameters(lin_solver::Tuple)
  return SSNParameters((
    maxnit=100, maxlit=50, tolsub=1e-7,
    ls_method=:strongwolfe, ls_tau=0.9, ls_c1=1e-4, ls_c2=0.9, stagit=2,
    lin_solver=lin_solver)...
  )
end
mutable struct ALMParameters{Tf<:AbstractFloat, Ti<:Integer}
  maxait  ::  Ti
  stoptol  ::  Tf
  maxn  ::  Ti
  maxtime  ::  Tf # max runtime
  stoptol_i1  ::  Tf # intermediate stop tol1
  stoptol_i2  ::  Tf # intermediate stop tol2
  BXratio :: Tf # ratio of X to B sigma
end
function ALMParameters()
  return ALMParameters((maxait = 100, stoptol = 1e-6, maxn = 50_000, maxtime=3600.0, stoptol_i1=1e-3, stoptol_i2=1e-6, BXratio=1.0)...)
end
