using SparseMatricesCSR

function default_options()
  options = Dict()
  # shared
  options[:maxtime] = 7_200
  options[:maxiter] = 100_000
  options[:max_oa_iter] = 1_000
  options[:kfactor] = 1
  options[:kfactor_scale] = 1
  options[:kfactor] = 3
  options[:kfactor_scale] = 1.0
  options[:kfactor_lowtol_iters] = 10
  options[:kfactor_lowtol] = 1e-3
  
  # osqp
  options[:scaled_termination] = false
  options[:eps_abs] = 1e-3 # lo-precision absolute: NOT gap (compl auto satisfied by osqp) #! main control
  options[:eps_rel] = 1e-5 # lo-precision relative: NOT gap (compl auto satisfied by osqp)
  options[:pinfeas_tol_lo] = 1e-3 # lo-precision #! preceded by eps_abs
  options[:dinfeas_tol_lo] = 1e-3 # lo-precision #! preceded by eps_abs
  options[:polish] = false
  options[:polish_delta] = 1e-7 # regularization parameter for polishing
  options[:polish_refine_iter] = 10 # number of polish steps
    # # overall 1e-7
    # options[:eps_rel] = 1e-10
    # options[:eps_abs] = 1e-10
    # options[:pinfeas_tol_lo] = 1e-3
    # options[:dinfeas_tol_lo] = 1e-3
    # # overall 1e-2
    # options[:eps_rel] = 1e-3
    # options[:eps_abs] = 1e-3
    # options[:pinfeas_tol_lo] = 1e-10
    # options[:dinfeas_tol_lo] = 1e-10
    # # overall 1e-2
    # options[:eps_rel] = 1e-5
    # options[:eps_abs] = 1e-3
    # options[:pinfeas_tol_lo] = 1e-100
    # options[:dinfeas_tol_lo] = 1e-100
  
  # gurobi
  options[:method] = 2 # barrier
  options[:pinfeas_tol_hi] = 1e-6 # hi-precision
  options[:dinfeas_tol_hi] = 1e-6 # hi-precision
  options[:relgap_tol_hi] = 1e-6 # hi-precision
  options[:crossover] = 0 #! no crossover (barrier onlly)
  options[:presolve] = -1 #! -1 default
  options[:scaleflag] = -1 # -1 default; 2 geometric mean
  options[:numericfocus] = 0 # 0 default; 3 most numerically conscious
  options[:nthreads] = 4 # BLAS.get_num_threads() # match p-ALM
  return options
end

function calcs_manual!(
  res::Dict, obj::FunctionValue, PI::ProblemInstance,
  x::Vector{Tf}, y::Vector{Tf}, z::Vector{Tf}, t::Vector{Tf},
  lambda::Vector{Tf}, mu::Vector{Tf}
) where Tf<:AbstractFloat
  m, n = PI.m, PI.n
  M = sum(m)
  L = length(m)
  
  # variables
  res[:x] = copy(x)
  res[:t] = copy(t)
  res[:y] = [zeros(m[l]) for l in 1:L]
  res[:z] = [zeros(m[l]) for l in 1:L]
  res[:lambda] = [zeros(m[l]) for l in 1:L]
  res[:mu] = copy(mu)
  for l in 1:L
    lidx = sum(m[1:l-1])+1:sum(m[1:l])
    res[:y][l] .= copy(y[lidx])
    res[:z][l] .= copy(z[lidx])
    res[:lambda][l] .= copy(lambda[lidx])
  end

  # primal infeasibility
  res[:maxksum] = [maxksum(res[:y][l], PI.k[l]) for l in 1:L]
  res[:maxksum_infeas] = maximum(max.(res[:maxksum] .- PI.r, 0.0))
  res[:maxksumrhs] = PI.r
  res[:cvar] = [maxksum(res[:y][l], PI.k[l])/PI.k[l] for l in 1:L]
  res[:cvar_rhs] = PI.r./PI.k
  res[:cvar_infeas] = maximum(max.(res[:cvar] .- res[:cvar_rhs], 0.0))    
  res[:yAxb_infeas] = maximum(max.([norm(res[:y][l] .- (PI.A[l] * res[:x] + PI.b[l]), Inf)/(1+norm(PI.b[l])) for l in 1:PI.L]))
  res[:xlohi_infeas] = max(norm(max.(res[:x] .- PI.xhi, 0.0), Inf)/(1+norm(res[:x], Inf)), norm(max.(PI.xlo .- res[:x], 0.0), Inf)/(1+norm(res[:x], Inf)))
  res[:pinfeas] = max(res[:yAxb_infeas], res[:xlohi_infeas], res[:cvar_infeas])
  
  # primal and dual objective
  Atlampmu = PI.A' * res[:lambda] .+ res[:mu]
  Btlam = [zeros(m[l]) for l in eachindex(m)]
  # for l in 1:PI.L
  #   if isa(PI.Bt[l], Function)
  #     @views PI.Bt[l](Btlam[l], res[:lambda][l])
  #   else
  #     @views Btlam[l] = PI.Bt[l] * res[:lambda][l]
  #   end
  # end
  pobj!(obj, res[:x], PI.C, PI.c, zeros(PI.n))
  Xunconstr = false
  tol = 1e-6
  dobj!(
    obj, PI.objtype,
    res[:lambda], res[:mu], zeros(PI.n), zeros(PI.n),
    Atlampmu, Btlam,
    PI.b, PI.r, PI.k,
    PI.xlo, PI.xhi, PI.xloidx, PI.xhiidx, PI.Xunconstr,
    PI.Cinv, PI.c,
    tol
  )
  res[:pobj] = obj.v
  res[:dobj] = obj.d
  res[:relgap] = (res[:pobj] - res[:dobj]) / (1.0 + abs(res[:pobj]))
  
  # dual infeasibility
  res[:dinfeas] = norm(obj.g .- Atlampmu, Inf)/(1+norm(obj.g, Inf))

  nothing
end

function solve_cvar_reduced_osqp(PI::ProblemInstance, options::Dict, destroyinput::Bool=true)
  """
  min   f(x)
  s.t.  sum(z)/m <= t*(1-alpha)    : c1
        z >= t1 + (Ax+b) - (r/k)1  : c2
        z >= 0                     : c3
        t >= 0                     : c4 -- allowed because of c1
  
  as

  min   0.5 x' * P * x + q' * x
  s.t.  l <= Ax <= u

  note: this formulation is slightly tighter than the `direct LP-existence formulation`
        since `r/k-s >= 0, s ∈ ℝ` would give only `s < 0 <= r/k`
  """
  start_time = time()

  #
  # initialize
  #

  # dims
  n = PI.n
  m = PI.m
  M = sum(PI.m)
  L = length(PI.m)
  nn = PI.n + M + L # X = [x, z, t]
  mm = L + M + M + PI.n + L # sum(z) <= t; t + y - r/k <= z; z >= 0; l <= x <= u; 0 <= t
  println("(nn, mm) = ($nn, $mm)")
  
  # alpha
  alpha = (1 .- Rational.(PI.k) ./ Rational.(PI.m))

  #
  # P matrix
  #

  if PI.objtype == 2
    P = [
      PI.C              spzeros(n, nn-n);
      spzeros(nn-n, n)  spzeros(nn-n, nn-n)
    ]
  else
    P = [
      spzeros(n,n)      spzeros(n, nn-n);
      spzeros(nn-n, n)  spzeros(nn-n, nn-n)
    ]
  end
  
  #
  # q vector
  #

  q = zeros(nn);
  q[1:PI.n] .= PI.c

  #
  # A matrix
  #

  aaa = hcat([[spzeros(sum(m[1:l-1])); ones(m[l]); spzeros(sum(m[l+1:end]))] for l in eachindex(m)]...)
  bbb = hcat([[spzeros(sum(m[1:l-1])); ones(m[l]) ./ m[l]; spzeros(sum(m[l+1:end]))] for l in eachindex(m)]...)'
  A = [
  # x                  z                      t
    spzeros(L, PI.n)   bbb                    Diagonal(-(1 .- alpha))   # c1: sum(zl)/m - tl * (1-alpha) <= 0
    vcat(PI.A...)      -spdiagm(0=>ones(M))   aaa;                      # c2: zl >= tl .* 1l + (Al*x + bl) - 1l .* (rl/kl)
    spzeros(M, PI.n)   +I(M)                  spzeros(M, L);            # c3: 0 <= zl <= Inf
    +I(PI.n)           spzeros(PI.n, M)       spzeros(PI.n, L);         # c4: l <= x <= u
    spzeros(L, PI.n)   spzeros(L, M)          +I(L);                    # c5: 0 <= t <= Inf
  ]
  println("A time = $(time()-start_time)")

  #
  # bounds l and u
  #

  # l
  l = zeros(mm)
  u = zeros(mm)
  
  # constraint 1: [L] sum(zl)/m - tl * (1-alpha) <= 0
  roffset = 0
  @views l[roffset+1:roffset+L] .= -Inf
  @views u[roffset+1:roffset+L] .= 0.0
  roffset += L
  
  # constraint 2: [M] Al*x - zl + tl .* 1l <= 1l .* (rl/kl) - bl
  @views l[roffset+1:roffset+M] .= -Inf
  @views u[roffset+1:roffset+M] .= vcat([PI.r[l] ./ PI.k[l] .- PI.b[l] for l in eachindex(m)]...)
  roffset += M

  # constraint 3: [M] 0 <= zl <= Inf
  @views l[roffset+1:roffset+M] .= 0.0
  @views u[roffset+1:roffset+M] .= +Inf
  roffset += M

  # constraint 4: l <= x <= u
  @views l[roffset+1:roffset+n] .= PI.xlo
  @views u[roffset+1:roffset+n] .= PI.xhi
  roffset += n

  # constraint 5: [L] 0 <= t <= Inf
  @views l[roffset+1:roffset+L] .= 0.0
  @views u[roffset+1:roffset+L] .= +Inf
  roffset += L

  # construct model
  model = OSQP.Model()
  @time OSQP.setup!(
    model; P=P, q=q, A=A, l=l, u=u,
    eps_prim_inf=options[:pinfeas_tol_lo],
    eps_dual_inf=options[:dinfeas_tol_lo],
    eps_rel=options[:eps_rel], # e_prim = e_abs + e_rel*max(|Ax|,|z|) and e_dual = e_abs + e_rel*max(|Px|,|Aty|,|q|)
    eps_abs=options[:eps_abs],
    time_limit=options[:maxtime],
    max_iter=options[:maxiter],
    # linsys_solver=options[:linsys_solver], #! needs PARDISO MKL
    linsys_solver="qdldl",
    scaled_termination=options[:scaled_termination],
    polish=options[:polish],
    delta=options[:polish_delta],
    polish_refine_iter=options[:polish_refine_iter],
  )
  nz = nnz(P) + nnz(A)
  A = nothing; P = nothing; aaa = nothing; bbb = nothing; GC.gc()
  if destroyinput
    PI = nothing; GC.gc()
  end
  println("removed A")
  println("model initialization = $(time()-start_time)")
  # solve model
  results = OSQP.solve!(model)
  
  # result
  xidx = 1:n
  zidx = n+1:n+M
  tidx = n+M+1:n+M+L
  lambdaidx = L+1:L+M
  muidx = 2M+L+1:2M+L+n
  res = Dict()
  res[:x] = results.x[xidx]
  res[:normx] = norm(res[:x])
  # res[:y] = vcat([PI.A[l] * res[:x] .+ PI.b[l] for l in 1:L]...)
  res[:z] = results.x[zidx]
  res[:t] = results.x[tidx]
  res[:lambda] = -results.y[lambdaidx]
  res[:mu] = -results.y[muidx]
  res[:walltime] = results.info.run_time # setup time included in OSQP Table 2 timing
  # res[:walltime] = results.info.solve_time + results.info.polish_time
  res[:status] = (results.info.status == :Solved ? 1 : 0)
  res[:iter] = results.info.iter
  res[:nnz] = nz
  model = nothing; GC.gc()
  return res
end

function solve_cvar_reduced_gurobi(PI::ProblemInstance, options::Dict, destroyinput::Bool=true, scaling::Symbol=:k1)
  """
  `:m` scaling
  min   f(x)
  s.t.  sum(z)/m <= t*(1-alpha)    : c1
        z >= t1 + (Ax+b) - (r/k)1  : c2
        z >= 0                     : c3
        t >= 0                     : c4 -- allowed because of c1
  note: this formulation is slightly tighter than the `direct LP-existence formulation`
        since `r/k-s >= 0, s ∈ ℝ` would give only `s < 0 <= r/k`
  `:k1` scaling (standard)
  min   f(x)
  s.t.  t + sum(z)/k <= 0          : c1
        z >= (Ax+b) - t1 - (r/k)1  : c2 <==> -b + (r/k)1 >= Ax - z - t1
        z >= 0                     : c3
  `:k2` scaling (standard)
  min   f(x)
  s.t.  t + sum(z)/k <= r/k        : c1
        z >= (Ax+b) - t1           : c2 <==> -b >= Ax - z - t1
        z >= 0                     : c3
  """
  start_time = time()
  
  #
  # initialize
  #

  # dims
  n, m, k = PI.n, PI.m, PI.k
  M = sum(PI.m)
  L = length(PI.m)
  nn = PI.n + M + L # X = [x, z, t]
  mm = M + L # t + (Ax+b) - r/k <= z; sum(z) <= t
  nn_Cint = Cint(nn)
  nn_Csize_t = Base.Csize_t(nn)
  xidx = 1:n
  zidx = n+1:n+M
  tidx = n+M+1:n+M+L
  
  # alpha
  if scaling == :m
    alpha = (1 .- Rational.(PI.k) ./ Rational.(PI.m))
  end

  #
  # A matrix
  #

  if scaling == :m
    aaa = hcat([[spzeros(sum(m[1:l-1])); ones(m[l]); spzeros(sum(m[l+1:end]))] for l in eachindex(m)]...)
    bbb = hcat([[spzeros(sum(m[1:l-1])); ones(m[l]) ./ m[l]; spzeros(sum(m[l+1:end]))] for l in eachindex(m)]...)'
    Acsr = [
    # x                  z                      t
      spzeros(L, PI.n)   bbb                    Diagonal(-(1 .- alpha))   # c1: sum(zl)/m - tl * (1-alpha) <= 0
      vcat(PI.A...)      -spdiagm(0=>ones(M))   aaa;                      # c2: zl >= tl .* 1l + (Al*x + bl) - 1l .* (rl/kl)
    ]
  elseif scaling == :k1
    bbb = hcat([[spzeros(sum(m[1:l-1])); ones(m[l]) ./ k[l]; spzeros(sum(m[l+1:end]))] for l in eachindex(m)]...)'
    aaa = hcat([[spzeros(sum(m[1:l-1])); -ones(m[l]); spzeros(sum(m[l+1:end]))] for l in eachindex(m)]...)
    Acsr = [
    # x                  z                      t
      spzeros(L, PI.n)   bbb                    I     # c1: sum(zl)/kl + tl <= 0
      vcat(PI.A...)      -spdiagm(0=>ones(M))   aaa;  # c2: (Al*x + bl) - rl/kl - zl - tl <= 0
    ]
  elseif scaling == :k2
    # note: same matrix as `:k1`
    bbb = hcat([[spzeros(sum(m[1:l-1])); ones(m[l]) ./ k[l]; spzeros(sum(m[l+1:end]))] for l in eachindex(m)]...)'
    aaa = hcat([[spzeros(sum(m[1:l-1])); -ones(m[l]); spzeros(sum(m[l+1:end]))] for l in eachindex(m)]...)
    Acsr = [
    # x                  z                      t
      spzeros(L, PI.n)   bbb                    I     # c1: sum(zl)/kl + tl - rl/kl <= 0
      vcat(PI.A...)      -spdiagm(0=>ones(M))   aaa;  # c2: (Al*x + bl) - zl - tl <= 0
    ]
  end
  c1idx = 1:L
  c2idx = L+1:L+M
  Acsr = SparseMatricesCSR.SparseMatrixCSR(transpose(sparse(transpose(Acsr))));
  @assert(Acsr.m == mm)
  # writedlm("A_reduced.csv", A)
  aaa = nothing; bbb = nothing; GC.gc()
  println("collected A constraint")

  #
  # environment
  #

  env_p = Ref{Ptr{Cvoid}}()
  logfile = C_NULL
  # logfile = "cvar.lp"
  error = GRBloadenv(env_p, logfile)
  @assert error == 0
  env = env_p[]
  error = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 1)
  @assert error == 0
  
  # parameters
  error = GRBsetdblparam(env, "TimeLimit", options[:maxtime])
  @assert error == 0
  method = (haskey(options, :method) ? options[:method] : 2) # 2 == barrier
  @assert error == 0
  error = GRBsetintparam(env, "Method", Cint(method))
  @assert error == 0
  error = GRBsetintparam(env, "BarIterLimit", options[:maxiter])
  @assert error == 0
  error = GRBsetdblparam(env, "FeasibilityTol", options[:pinfeas_tol_hi])
  @assert error == 0
  error = GRBsetdblparam(env, "OptimalityTol", options[:dinfeas_tol_hi])
  @assert error == 0
  error = GRBsetdblparam(env, "BarConvTol", options[:relgap_tol_hi])
  @assert error == 0
  error = GRBsetintparam(env, "Threads", options[:nthreads])
  @assert error == 0
  error = GRBsetintparam(env, "Presolve", options[:presolve])
  @assert error == 0
  error = GRBsetintparam(env, "Crossover", options[:crossover])
  @assert error == 0
  error = GRBsetintparam(env, "ScaleFlag", options[:scaleflag]) # https://www.gurobi.com/documentation/10.0/refman/scaleflag.html#parameter:ScaleFlag
  @assert error == 0
  
  #
  # model
  #

  model_p = Ref{Ptr{Cvoid}}()
  modelname = "cvar"
  error = GRBnewmodel(env, model_p, modelname, 0, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)
  @assert error == 0
  model = model_p[]
  error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE)
  @assert error == 0

  #
  # bounds
  #

  xzt_lb = zeros(Cdouble, nn)
  xzt_ub = zeros(Cdouble, nn)
  
  # xlo <= x <= xhi
  xzt_lb[xidx] .= PI.xlo
  xzt_ub[xidx] .= PI.xhi

  # 0 <= z
  xzt_lb[zidx] .= 0.0
  xzt_ub[zidx] .= +Inf

  # 0 <= t
  if scaling == :m
    # xzt_lb[tidx] .= 0.0
    # xzt_ub[tidx] .= +Inf
  elseif scaling == :k1 || scaling == :k2
    xzt_lb[tidx] .= -Inf
    xzt_ub[tidx] .= +Inf
  end

  #
  # linear objective
  #

  xzt_lobj = zeros(Cdouble, nn)
  xzt_lobj[xidx] .= PI.c
  
  #
  # variables
  #

  error = GRBXaddvars(
      model,      # : model
      nn_Csize_t, # : numvars
      0,          # : numnz
      C_NULL,     # : *vbeg
      C_NULL,     # : *vind
      C_NULL,     # : *vval
      xzt_lobj,   # : *obj
      xzt_lb,     # : *lb
      xzt_ub,     # : *ub
      C_NULL,     # : *vtype
      C_NULL      # : **varnames
  )
  @assert error == 0

  #
  # quadratic objective
  #

  if PI.objtype == 2
    qrow, qcol, qval = findnz(sparse(PI.C))
    error = GRBaddqpterms(
      model,          # GRBmodel : *model
      Cint(n),        # int      : numqnz
      Cint.(qrow.-1), # int      : *qrow
      Cint.(qcol.-1), # int      : *qcol
      0.5qval,        # double   : *qval #! why 0.5???
    )
    @assert error == 0
  end

  #
  # constraints
  #

  consense = fill(Cchar(0), Acsr.m);
  consense[c1idx] .= Cchar(GRB_LESS_EQUAL);
  consense[c2idx] .= Cchar(GRB_LESS_EQUAL);
  conname = fill(C_NULL, Acsr.m);
  if scaling == :m
    rhs = [
      zeros(Cdouble, L); # c1
      Cdouble.(vcat([PI.r[l]/PI.k[l] .* ones(m[l]) .- PI.b[l] for l in eachindex(m)]...)); # c2
    ]
  elseif scaling == :k1
    # k1
    rhs = [
      zeros(Cdouble, L); # c1
      Cdouble.(vcat([PI.r[l]./PI.k[l] .* ones(m[l]) .- PI.b[l] for l in eachindex(m)]...)); # c2
    ]
  elseif scaling == :k2
    # k2
    rhs = [
      PI.r./PI.k; # c1
      Cdouble.(vcat([-PI.b[l] for l in eachindex(m)]...)); # c2
    ]
  end
  error = GRBXaddconstrs(
    model,                         # GRBmodel   : *model
    Cint(Acsr.m),                  # int        : numconstrs
    Base.Csize_t(nnz(Acsr)),       # size_t     : numnz
    Base.Csize_t.(Acsr.rowptr.-1), # size_t     : *cbeg
    Cint.(Acsr.colval.-1),         # int        : *cind
    Cdouble.(Acsr.nzval),          # double     : *cval
    consense,                      # char       : sense
    rhs,                           # double     : *rhs
    conname,                       # const char : **constrname
  )
  @assert error == 0
  Acsr = nothing; GC.gc()
  println("collected Acsr constraint")
  
  #
  # update solver parameters
  #

  # update
  error = GRBupdatemodel(model)
  @assert error == 0
  println("model initialization = $(time()-start_time)")

  # # force model updates
  # GRBsetdblparam(env, "TimeLimit", 0.0)
  # solve_time = @elapsed error = GRBoptimize(model)
  # println("model setup = $(solve_time)")
  
  if destroyinput
    PI = nothing; GC.gc()
    println("removed PI")
  end

  #
  # solve
  #

  try
    solve_time = @elapsed error = GRBoptimize(model)
  catch e
    println(e)
    err = GRBterminate(model)
    println(err)
    err = GRBfreemodel(model)
    println(err)
    GRBfreeenv(env)
    println("freed model")
    GC.gc()
    return
  end
  println("model solve = $(solve_time)")
  
  #
  # recover solution
  #

  res = Dict();
  NumVars = Ref{Cint}();
  NumConstrs = Ref{Cint}();
  IterCount = Ref{Cdouble}(); # simplex iters
  BarIterCount = Ref{Cint}(); # barrier iters
  ObjVal = Ref{Cdouble}(); # primal obj
  ObjBound = Ref{Cdouble}(); # dual obj
  RunTime = Ref{Cdouble}(); # solver wall-clock time
  xzt = zeros(Cdouble, nn); # xzt
  nu_1 = zeros(Cdouble, L); # c1
  nu_2 = zeros(Cdouble, M); # c2
  mu = zeros(Cdouble, n); # xlo <= x <= xhi
  Status = Ref{Cint}(); # solved status
  DNumNZs = Ref{Cdouble}() # number of nonzero coeffs in linear constraints
  
  error = GRBgetintattr(model, "NumVars", NumVars);
  @assert error == 0
  error = GRBgetintattr(model, "NumConstrs", NumConstrs);
  @assert error == 0
  error = GRBgetintattr(model, "Status", Status);
  @assert error == 0
  if method == 2
    error = GRBgetintattr(model, "BarIterCount", BarIterCount);
    @assert error == 0
  else
    error = GRBgetdblattr(model, "IterCount", IterCount);
    @assert error == 0
  end
  error = GRBgetdblattr(model, "ObjVal", ObjVal);
  @assert error == 0
  # error = GRBgetdblattr(model, "ObjBound", ObjBound); #! barrier doesn't have this...
  # @assert error == 0
  error = GRBgetdblattr(model, "RunTime", RunTime);
  @assert error == 0
  GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, nn, xzt);
  @assert error == 0
  GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, 0, L, nu_1);
  @assert error == 0
  GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, L, M, nu_2);
  @assert error == 0
  GRBgetdblattrarray(model, GRB_DBL_ATTR_RC, 0, n, mu);
  @assert error == 0
  error = GRBgetdblattr(model, "DNumNZs", DNumNZs);
  @assert error == 0
  
  # populate result
  res[:nnz] = DNumNZs[]
  res[:numvar] = NumVars[]
  res[:numcon] = NumConstrs[]
  res[:walltime] = RunTime[]
  if method == 2
    res[:iter] = BarIterCount[]
  else
    res[:iter] = IterCount[]
  end
  res[:pobj] = ObjVal[]
  # res[:dobj] = ObjBound[]
  res[:x] = xzt[xidx]
  res[:normx] = norm(res[:x])
  res[:z] = xzt[zidx]
  res[:t] = xzt[tidx]
  # res[:nu_1] = -nu_1
  # res[:nu_2] = -nu_2
  res[:mu] = mu
  res[:lambda] = nu_2
  res[:status] = (Status[] == 2 ? 1 : 0)
  # res[:y] = vcat([PI.A[l] * res[:x] .+ PI.b[l] for l in 1:L]...)
  
  #
  # free model
  #
  
  # error = GRBwrite(model, "cvar_direct.mps");
  @assert error == 0
  error = GRBfreemodel(model)
  @assert error == 0
  GRBfreeenv(env)
  return res
end





#==========
 BENCHMARK
==========#
function get_partition(y::Vector, k::Union{Integer,Vector{Ti}}) where Ti<:Integer
  if eltype(y) <: Vector
    partition = Array{Array{Int64}}[]
    for i in eachindex(y)
      sig = sortperm(y[i], rev=true)
      k0bar, k1bar = get_k0k1(y[i][sig], k[i])
      alpha = sig[1:k0bar]
      beta = sig[k0bar+1:k1bar]
      gamma = sig[k1bar+1:end]
      push!(partition, [alpha, beta, gamma])
    end
  else
    sig = sortperm(y, rev=true)
    k0bar, k1bar = get_k0k1(y[sig], k)
    alpha = sig[1:k0bar]
    beta = sig[k0bar+1:k1bar]
    gamma = sig[k1bar+1:end]
    partition = [alpha, beta, gamma]
  end
  return partition
end
function solve_cvar_reduced_gurobi_benchmark(PI::ProblemInstance, partition::Vector, options::Dict, destroyinput::Bool=true)
  """
  min   f(x)
  s.t.  sum( t1_alpha + (Ax+b)_alpha - (r/k)1_alpha )/m <= t*(1-alpha) : c1
        0 <= t1_alpha + (Ax+b)_alpha - (r/k)1_alpha                    : c2
        0 == t1_beta + (Ax+b)_beta - (r/k)1_beta                       : c3
        0 >= t1_gamma + (Ax+b)_gamma - (r/k)1_gamma                    : c4
        t >= 0                                                         : c5 -- allowed because of c1
  
  note: this formulation is slightly tighter than the `direct LP-existence formulation`
        since `r/k-s >= 0, s ∈ ℝ` would give only `s < 0 <= r/k`; it sets z_betagamma == 0
  """
  start_time = time()
  
  #
  # initialize
  #

  # dims
  n, m = PI.n, PI.m
  M = sum(PI.m)
  L = length(PI.m)
  
  nn = PI.n + L # X = [x, zalpha, t]
  mm = L + M # sum(zalpha) <= t; t + (Ax+b) - r/k >=/==/<= 0
  nn_Cint = Cint(nn)
  nn_Csize_t = Base.Csize_t(nn)
  xidx = 1:n
  tidx = n+1:n+L
  
  # alpha
  alpha = (1 .- Rational.(PI.k) ./ Rational.(PI.m))

  # alpha partition
  AP = [partition[l][1] for l in 1:L] # alpha partition
  lAP = [length(ap) for ap in AP] # length of alpha partition components
  nAP = sum(lAP) # size of alpha partition
  
  # beta partition
  BP = [partition[l][2] for l in 1:L]
  lBP = [length(bp) for bp in BP]
  nBP = sum(length(bp) for bp in BP)
  
  # gamma partition
  GP = [partition[l][3] for l in 1:L]
  lGP = [length(gp) for gp in GP]
  nGP = sum(length(gp) for gp in GP)
  
  #
  # A matrix
  #
  Aalpha = vcat([PI.A[l][AP[l],:] for l in 1:L]...)
  Abeta = vcat([PI.A[l][BP[l],:] for l in 1:L]...)
  Agamma = vcat([PI.A[l][GP[l],:] for l in 1:L]...)
  sumAalpha = vcat([(ones(lAP[l]) ./ m[l])' * PI.A[l][AP[l],:] for l in 1:L]...)
  # sumAalpha = vcat([(ones(lAP[l]))' * PI.A[l][AP[l],:] for l in 1:L]...)
  indalpha = hcat([[spzeros(sum(lAP[1:l-1])); ones(lAP[l]); spzeros(sum(lAP[l+1:end]))] for l in 1:L]...)
  indbeta = hcat([[spzeros(sum(lBP[1:l-1])); ones(lBP[l]); spzeros(sum(lBP[l+1:end]))] for l in 1:L]...)
  indgamma = hcat([[spzeros(sum(lGP[1:l-1])); ones(lGP[l]); spzeros(sum(lGP[l+1:end]))] for l in 1:L]...)
  
  # x           # t                          # description
  Acsr = [
    sumAalpha   spdiagm((lAP .- PI.k) ./ m)  # c1: sum ( (Ax+b)_alpha + t - r/k ) / m <= t*(1-tau) where tau="alpha_threshold"
    # sumAalpha   spdiagm((lAP .- PI.k))       # c1: sum ( (Ax+b)_alpha + t - r/k ) / m <= t*(1-tau) where tau="alpha_threshold"
    Aalpha      indalpha                     # c2: (Ax+b)_alpha + t - r/k >= 0
    Abeta       indbeta                      # c3: (Ax+b)_alpha + t - r/k == 0
    Agamma      indgamma                     # c4: (Ax+b)_alpha + t - r/k <= 0
  ]
  sumAalpha = nothing; Aalpha = nothing; Abeta = nothing; Agamma = nothing;
  indalpha = nothing; indbeta = nothing; indgamma = nothing;
  GC.gc()
  c1idx = 1:L
  c2idx = L+1:L+nAP
  c3idx = L+nAP+1:L+nAP+nBP
  c4idx = L+nAP+nBP+1:L+nAP+nBP+nGP
  Acsr = SparseMatricesCSR.SparseMatrixCSR(transpose(sparse(transpose(Acsr))));
  @assert(Acsr.m == mm)
  # writedlm("A_reduced.csv", A)
  println("collected A constraint")

  #
  # environment
  #

  env_p = Ref{Ptr{Cvoid}}()
  logfile = C_NULL
  # logfile = "cvar.lp"
  error = GRBloadenv(env_p, logfile)
  @assert error == 0
  env = env_p[]
  error = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 1)
  @assert error == 0
  
  # parameters
  error = GRBsetdblparam(env, "TimeLimit", options[:maxtime])
  @assert error == 0
  method = (haskey(options, :method) ? options[:method] : 2) # 2 == barrier
  @assert error == 0
  error = GRBsetintparam(env, "Method", Cint(method))
  @assert error == 0
  error = GRBsetintparam(env, "BarIterLimit", options[:maxiter])
  @assert error == 0
  error = GRBsetdblparam(env, "FeasibilityTol", options[:pinfeas_tol_hi])
  @assert error == 0
  error = GRBsetdblparam(env, "OptimalityTol", options[:dinfeas_tol_hi])
  @assert error == 0
  error = GRBsetdblparam(env, "BarConvTol", options[:relgap_tol_hi])
  @assert error == 0
  error = GRBsetintparam(env, "Threads", options[:nthreads])
  @assert error == 0
  error = GRBsetintparam(env, "Presolve", options[:presolve])
  @assert error == 0
  error = GRBsetintparam(env, "Crossover", options[:crossover])
  @assert error == 0
  error = GRBsetintparam(env, "ScaleFlag", options[:scaleflag]) # https://www.gurobi.com/documentation/10.0/refman/scaleflag.html#parameter:ScaleFlag
  @assert error == 0
  error = GRBsetintparam(env, "NumericFocus", options[:numericfocus]) # https://www.gurobi.com/documentation/current/refman/numericfocus.html#parameter:NumericFocus
  @assert error == 0
  
  #
  # model
  #

  model_p = Ref{Ptr{Cvoid}}()
  modelname = "cvar"
  error = GRBnewmodel(env, model_p, modelname, 0, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)
  @assert error == 0
  model = model_p[]
  error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE)
  @assert error == 0

  #
  # bounds
  #

  xzt_lb = zeros(Cdouble, nn)
  xzt_ub = zeros(Cdouble, nn)
  
  # xlo <= x <= xhi
  xzt_lb[xidx] .= PI.xlo
  xzt_ub[xidx] .= PI.xhi

  # 0 <= t
  xzt_lb[tidx] .= 0.0
  xzt_ub[tidx] .= +Inf

  #
  # linear objective
  #

  xzt_lobj = zeros(Cdouble, nn)
  xzt_lobj[xidx] .= PI.c
  
  #
  # variables
  #

  error = GRBXaddvars(
      model,      # : model
      nn_Csize_t, # : numvars
      0,          # : numnz
      C_NULL,     # : *vbeg
      C_NULL,     # : *vind
      C_NULL,     # : *vval
      xzt_lobj,   # : *obj
      xzt_lb,     # : *lb
      xzt_ub,     # : *ub
      C_NULL,     # : *vtype
      C_NULL      # : **varnames
  )
  @assert error == 0

  #
  # quadratic objective
  #

  if PI.objtype == 2
    qrow, qcol, qval = findnz(sparse(PI.C))
    error = GRBaddqpterms(
      model,          # GRBmodel : *model
      Cint(n),        # int      : numqnz
      Cint.(qrow.-1), # int      : *qrow
      Cint.(qcol.-1), # int      : *qcol
      0.5qval,        # double   : *qval #! why 0.5???
    )
    @assert error == 0
  end

  #
  # constraints
  #

  consense = fill(Cchar(0), Acsr.m)
  consense[c1idx] .= Cchar(GRB_LESS_EQUAL)
  consense[c2idx] .= Cchar(GRB_GREATER_EQUAL)
  consense[c3idx] .= Cchar(GRB_EQUAL)
  consense[c4idx] .= Cchar(GRB_LESS_EQUAL)
  conname = fill(C_NULL, Acsr.m)
  rhs = [
    Cdouble.(((lAP .* PI.r ./ PI.k) .- [sum(PI.b[l][AP[l]]) for l in 1:L]) ./ m);                     # c1
    # Cdouble.(((lAP .* PI.r ./ PI.k) .- [sum(PI.b[l][AP[l]]) for l in 1:L]));                          # c1
    Cdouble.(vcat([(PI.r[l] / PI.k[l]) .* ones(Cdouble, lAP[l]) .- PI.b[l][AP[l]] for l in 1:L]...)); # c2
    Cdouble.(vcat([(PI.r[l] / PI.k[l]) .* ones(Cdouble, lBP[l]) .- PI.b[l][BP[l]] for l in 1:L]...)); # c3
    Cdouble.(vcat([(PI.r[l] / PI.k[l]) .* ones(Cdouble, lGP[l]) .- PI.b[l][GP[l]] for l in 1:L]...)); # c4
  ]
  error = GRBXaddconstrs(
    model,                         # GRBmodel   : *model
    Cint(Acsr.m),                  # int        : numconstrs
    Base.Csize_t(nnz(Acsr)),       # size_t     : numnz
    Base.Csize_t.(Acsr.rowptr.-1), # size_t     : *cbeg
    Cint.(Acsr.colval.-1),         # int        : *cind
    Cdouble.(Acsr.nzval),          # double     : *cval
    consense,                      # char       : sense
    rhs,                           # double     : *rhs
    conname,                       # const char : **constrname
  )
  @assert error == 0
  Acsr = nothing; GC.gc()
  println("collected Acsr constraint")
  
  #
  # update solver parameters
  #

  # update
  error = GRBupdatemodel(model)
  @assert error == 0
  println("model initialization = $(time()-start_time)")

  # # force model updates
  # GRBsetdblparam(env, "TimeLimit", 0.0)
  # solve_time = @elapsed error = GRBoptimize(model)
  # println("model setup = $(solve_time)")
  
  if destroyinput
    PI = nothing; GC.gc()
    println("removed PI")
  end

  #
  # solve
  #

  try
    solve_time = @elapsed error = GRBoptimize(model)
  catch e
    println(e)
    err = GRBterminate(model)
    println(err)
    err = GRBfreemodel(model)
    println(err)
    GRBfreeenv(env)
    println("freed model")
    GC.gc()
    return
  end
  println("model solve = $(solve_time)")
  
  #
  # recover solution
  #

  res = Dict();
  NumVars = Ref{Cint}();
  NumConstrs = Ref{Cint}();
  IterCount = Ref{Cdouble}(); # simplex iters
  BarIterCount = Ref{Cint}(); # barrier iters
  ObjVal = Ref{Cdouble}(); # primal obj
  ObjBound = Ref{Cdouble}(); # dual obj
  RunTime = Ref{Cdouble}(); # solver wall-clock time
  xzt = zeros(Cdouble, nn); # xzt
  # nu_1 = zeros(Cdouble, L); # c1
  # nu_2 = zeros(Cdouble, M); # c2
  # mu = zeros(Cdouble, n); # xlo <= x <= xhi
  Status = Ref{Cint}(); # solved status
  DNumNZs = Ref{Cdouble}() # number of nonzero coeffs in linear constraints
  
  error = GRBgetintattr(model, "NumVars", NumVars);
  @assert error == 0
  error = GRBgetintattr(model, "NumConstrs", NumConstrs);
  @assert error == 0
  error = GRBgetintattr(model, "Status", Status);
  @assert error == 0
  if method == 2
    error = GRBgetintattr(model, "BarIterCount", BarIterCount);
    @assert error == 0
  else
    error = GRBgetdblattr(model, "IterCount", IterCount);
    @assert error == 0
  end
  error = GRBgetdblattr(model, "ObjVal", ObjVal);
  @assert error == 0
  # error = GRBgetdblattr(model, "ObjBound", ObjBound); #! barrier doesn't have this...
  # @assert error == 0
  error = GRBgetdblattr(model, "RunTime", RunTime);
  @assert error == 0
  GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, nn, xzt);
  @assert error == 0
  # GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, 0, L, nu_1);
  # @assert error == 0
  # GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, L, M, nu_2);
  # @assert error == 0
  # GRBgetdblattrarray(model, GRB_DBL_ATTR_RC, 0, n, mu);
  @assert error == 0
  error = GRBgetdblattr(model, "DNumNZs", DNumNZs);
  @assert error == 0
  
  # populate result
  res[:nnz] = DNumNZs[]
  res[:numvar] = NumVars[]
  res[:numcon] = NumConstrs[]
  res[:walltime] = RunTime[]
  if method == 2
    res[:iter] = BarIterCount[]
  else
    res[:iter] = IterCount[]
  end
  res[:pobj] = ObjVal[]
  # res[:dobj] = ObjBound[]
  res[:x] = xzt[xidx]
  res[:t] = xzt[tidx]
  # res[:nu_1] = -nu_1
  # res[:nu_2] = -nu_2
  # res[:mu] = mu
  # res[:lambda] = nu_2
  res[:status] = (Status[] == 2 ? 1 : 0)
  # res[:y] = vcat([PI.A[l] * res[:x] .+ PI.b[l] for l in 1:L]...)
  
  #
  # free model
  #
  
  # error = GRBwrite(model, "cvar_direct.mps");
  @assert error == 0
  error = GRBfreemodel(model)
  @assert error == 0
  GRBfreeenv(env)
  return res
end
