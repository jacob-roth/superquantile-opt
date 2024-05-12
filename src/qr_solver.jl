function solve_qr_gurobi(PI::ProblemInstance, options::Dict)
  solve_qr_gurobi(PI.A[1], PI.b[1], PI.k[1], options)
end
function solve_qr_gurobi(X::AbstractMatrix, y::AbstractVector, k::Integer, options::Dict)
  """
  tau := (m-k) / m
  min  1/m sum( (1-tau) * zneg  + tau * zpos)
  s.t. zpos - zneg = y - X * theta
       zpos >= 0
       zneg >= 0
  """
  start_time = time()
  
  #
  # initialize
  #

  # dims
  m, n = size(X)
  nn = n + m + m # [x, zpos, zneg]
  mm = m # zpos - zneg = y - X * x; x == theta
  xidx = 1:n
  zposidx = n+1:n+m
  znegidx = n+m+1:n+2m
  
  # alpha
  tau = (m-k) / m

  #
  # A matrix
  #

  Acsr = [
  # x   zpos                   zneg
    X   +spdiagm(0=>ones(m))   -spdiagm(0=>ones(m))
  ]
  Acsr = SparseMatricesCSR.SparseMatrixCSR(transpose(sparse(transpose(Acsr))));
  @assert(Acsr.m == mm)
  
  #
  # environment
  #

  env_p = Ref{Ptr{Cvoid}}()
  logfile = C_NULL
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

  xzz_lb = zeros(Cdouble, nn)
  xzz_ub = zeros(Cdouble, nn)
  
  # xlo <= x <= xhi
  xzz_lb[xidx] .= -Inf
  xzz_ub[xidx] .= +Inf

  # 0 <= zpos
  xzz_lb[zposidx] .= 0.0
  xzz_ub[zposidx] .= +Inf

  # 0 <= zneg
  xzz_lb[znegidx] .= 0.0
  xzz_ub[znegidx] .= +Inf

  #
  # linear objective
  #

  xzz_lobj = zeros(Cdouble, nn)
  xzz_lobj[zposidx] .= tau
  xzz_lobj[znegidx] .= 1-tau
  
  #
  # variables
  #

  error = GRBXaddvars(
      model,      # : model
      Base.Csize_t(nn), # : numvars
      0,          # : numnz
      C_NULL,     # : *vbeg
      C_NULL,     # : *vind
      C_NULL,     # : *vval
      xzz_lobj,   # : *obj
      xzz_lb,     # : *lb
      xzz_ub,     # : *ub
      C_NULL,     # : *vtype
      C_NULL      # : **varnames
  )
  @assert error == 0

  #
  # constraints
  #

  consense = fill(Cchar(0), Acsr.m);
  consense .= Cchar(GRB_EQUAL);
  conname = fill(C_NULL, Acsr.m);
  rhs = y; # pointer
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
  MaxVio = Ref{Cdouble}(); # max of all violation
  ComplVio = Ref{Cdouble}(); # max compl violation
  ConstrVio = Ref{Cdouble}(); # max constr violation
  DualVio = Ref{Cdouble}(); # max dual violation
  xzz = zeros(Cdouble, nn); # xzz
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
  error = GRBgetdblattr(model, "MaxVio", MaxVio);
  @assert error == 0
  error = GRBgetdblattr(model, "ComplVio", ComplVio);
  @assert error == 0
  error = GRBgetdblattr(model, "ConstrVio", ConstrVio);
  @assert error == 0
  error = GRBgetdblattr(model, "DualVio", DualVio);
  @assert error == 0
  error = GRBgetdblattr(model, "ObjVal", ObjVal);
  @assert error == 0
  # error = GRBgetdblattr(model, "ObjBound", ObjBound); #! barrier doesn't have this...
  # @assert error == 0
  error = GRBgetdblattr(model, "RunTime", RunTime);
  @assert error == 0
  GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, nn, xzz);
  @assert error == 0
  # GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, 0, L, nu_1);
  # @assert error == 0
  # GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, L, M, nu_2);
  # @assert error == 0
  # GRBgetdblattrarray(model, GRB_DBL_ATTR_RC, 0, n, mu);
  # @assert error == 0
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
  res[:x] = xzz[xidx]
  res[:normx] = norm(res[:x])
  res[:zpos] = xzz[zposidx]
  res[:zneg] = xzz[znegidx]
  res[:status] = (Status[] == 2 ? 1 : 0)
  res[:pinfeas] = ConstrVio[]
  res[:dinfeas] = DualVio[]
  res[:relgap] = ComplVio[]
  
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
