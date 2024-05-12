
function solve_cvar_oa_reduced_gurobi(PI::ProblemInstance, options::Dict)
  """
  min   f(x)
  s.t.  t + sum(z_W)/k <= 0                : c1
        z_W >= (Ax+b)_W - (r/k)1_W - t1_W  : c2 <==> -b_W + (r/k)1_W >= (Ax)_W - z_W - t1_W
        z_W >= 0                           : c3
        where W ⊆ {1,...,m}
  """
  #
  # initialize
  #

  intermediate_kfactor = copy(options[:kfactor])
  intermediate_kfactor_scale = copy(options[:kfactor_scale])

  # output
  start_time = time(); # total walltime
  RunTime = Ref{Cdouble}(); # subproblem walltime
  IterCount = Ref{Cdouble}(); # simplex iters
  BarIterCount = Ref{Cint}(); # barrier iters
  res = Dict(); # output
  res[:times] = zeros(options[:max_oa_iter])
  res[:iters] = zeros(options[:max_oa_iter])
  retcode = 1

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
  alpha = (1 .- Rational.(PI.k) ./ Rational.(PI.m))

  #
  # A matrix
  #

  bbb = hcat([[spzeros(sum(m[1:l-1])); ones(m[l]) ./ k[l]; spzeros(sum(m[l+1:end]))] for l in eachindex(m)]...)'
  Acsr_c1 = [
  # x                  z                      t
    spzeros(L, PI.n)   bbb                    Diagonal(+ones(L))   # c1: sum(zl)/kl + tl <= 0
  ]
  c1idx = 1:L
  Acsr_c1 = SparseMatricesCSR.SparseMatrixCSR(transpose(sparse(transpose(Acsr_c1))));
  @assert(Acsr_c1.m == L)
  bbb = nothing; GC.gc()
  
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

  xzt_lb = zeros(Cdouble, nn)
  xzt_ub = zeros(Cdouble, nn)
  
  # xlo <= x <= xhi
  xzt_lb[xidx] .= PI.xlo
  xzt_ub[xidx] .= PI.xhi

  # 0 <= z
  xzt_lb[zidx] .= 0.0
  xzt_ub[zidx] .= 0.0

  # t free
  xzt_lb[tidx] .= -Inf
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

  consense_c1 = fill(Cchar(0), Acsr_c1.m)
  consense_c1 .= GRB_LESS_EQUAL
  conname_c1 = fill(C_NULL, Acsr_c1.m)
  rhs_c1 = zeros(Cdouble, L); # c1: t + sum(z_W)/k <= 0
  rhs_c2 = [Cdouble.(PI.r[l]/PI.k[l] .* ones(m[l]) .- PI.b[l]) for l in 1:L]; # c2: -b_W + (r/k)1_W >= (Ax)_W - z_W - t1_W <==> z_W >= (Ax+b)_W - (r/k)1_W - t1_W
  
  # add c1; save c2 for OA
  error = GRBXaddconstrs(
    model,                            # GRBmodel   : *model
    Cint(Acsr_c1.m),                  # int        : numconstrs
    Base.Csize_t(nnz(Acsr_c1)),       # size_t     : numnz
    Base.Csize_t.(Acsr_c1.rowptr.-1), # size_t     : *cbeg
    Cint.(Acsr_c1.colval.-1),         # int        : *cind
    Cdouble.(Acsr_c1.nzval),          # double     : *cval
    consense_c1,                      # char       : sense
    rhs_c1,                           # double     : *rhs
    conname_c1,                       # const char : **constrname
  )
  @assert error == 0

  # free memory
  Acsr_c1 = nothing; GC.gc()
  println("collected Acsr_c1 constraint")
  
  #
  # update solver parameters
  #

  # update
  error = GRBupdatemodel(model)
  @assert error == 0
  println("model initialization = $(time()-start_time)")

  #
  # prepare inactive
  #

  inactive = [ones(Bool, ml) for ml in m]
  active = [zeros(Bool, ml) for ml in m]
  active_prev = [zeros(Bool, ml) for ml in m]
  xzt = zeros(Cdouble, nn); # xzt
  vx = zeros(n)
  if PI.objtype == 1
    # for i in 1:n
    #   if PI.c[i] > 0
    #     vx[i] = isinf(PI.xlo[i]) ? 0.0 : PI.xlo[i]
    #   elseif PI.c[i] < 0
    #     vx[i] = isinf(PI.xhi[i]) ? 0.0 : PI.xhi[i]
    #   end
    # end
    vx .= 0.0
  elseif PI.objtype == 2
    vx .= -PI.Cinv * PI.c
  end
  vz = [zeros(ml) for ml in m]
  vt = zeros(L)
  gx = [zeros(ml) for ml in m]
  top_kfactor_idx = [1:Int(ceil(min(options[:kfactor] * PI.k[l], PI.m[l]))) for l in 1:L]
  sigma = [zeros(Int64, ml) for ml in m]
  oa_iter = 0
  
  #
  # solve
  #

  condict = Dict()
  zdict = Dict()
  # for l in 1:L
  #   condict[l] = Dict()
  # end
  order_added = 0
  tols = exp.(collect(range(log.(options[:kfactor_lowtol]), stop=log.(options[:relgap_tol_hi]), length=options[:kfactor_lowtol_iters])))
  while true

    # increment iteration
    oa_iter += 1
    
    # modify subproblem tolerance
    if oa_iter <= length(tols)
      tolsub = tols[oa_iter]
    else
      tolsub = min(options[:relgap_tol_hi], options[:pinfeas_tol_hi], options[:dinfeas_tol_hi])
    end
    # println(tolsub)
    error = GRBupdatemodel(model)
    @assert error == 0
    error = GRBsetdblparam(GRBgetenv(model), "FeasibilityTol", tolsub);
    # println(error)
    @assert error == 0
    error = GRBsetdblparam(GRBgetenv(model), "OptimalityTol", tolsub);
    # println(error)
    @assert error == 0
    error = GRBsetdblparam(GRBgetenv(model), "BarConvTol", tolsub);
    # println(error)
    @assert error == 0
    error = GRBupdatemodel(model)
    @assert error == 0
    mytol = Ref{Cdouble}(NaN)
    error = GRBgetdblparam(GRBgetenv(model), "BarConvTol", mytol)
    @assert error == 0

    # display subproblem info
    println("===========================================================\nsubproblem = $oa_iter, |active| = $([sum(active[l]) for l in 1:L]), subtol=$(mytol[])\n===========================================================\n")
    # println("
    # ==================================================================================================
    # subproblem = $oa_iter
    # |active| = $([sum(active[l]) for l in 1:L]), |active_prev| = $([sum(active_prev[l]) for l in 1:L])
    # |inactive| = $([sum(inactive[l]) for l in 1:L])
    # active = $([findall(active[l]) for l in 1:L]), active_prev = $([findall(active_prev[l]) for l in 1:L])
    # ==================================================================================================
    # ")

    # update active
    for l in 1:L
      if oa_iter > 1
        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, nn, xzt);
        @assert error == 0
        vx .= xzt[xidx]
        zlidx = zidx[1+sum(m[1:l-1]):sum(m[1:l])]
        vz[l] .= xzt[zlidx]
        vt[l] = xzt[tidx[l]]  
      end
      gx[l] .= PI.A[l] * vx .+ PI.b[l] .- PI.r[l]/PI.k[l] .- vz[l] .- vt[l]
      
      # union `k` largest
      intermediate_kfactor *= intermediate_kfactor_scale
      top_kfactor_idx = [1:Int(ceil(min(intermediate_kfactor * PI.k[l], PI.m[l]))) for l in 1:L]
      sigma[l] .= sortperm(gx[l], rev=true, alg=PartialQuickSort(top_kfactor_idx[l]))
      active[l][sigma[l][top_kfactor_idx[l]]] .= true
      #=
      # add min(remaining inactive, kfactor * k) elements
      sigma[l][inactive[l]] .= sortperm(gx[l][inactive[l]], rev=true, alg=PartialQuickSort(top_kfactor_idx[l]))
      if sum(inactive[l]) > 0
        top_kfactor_idx[l] = 1:min(length(active[l][sigma[l][inactive[l]]]), top_kfactor_idx[l][end])
        active[l][sigma[l][inactive[l]][top_kfactor_idx[l]]] .= true
      end
      =#
      inactive[l] .= .!active[l]
    end

    # set bounds on x for first iteration
    # if oa_iter == 1
    #   for i in 1:n
    #     @assert error == 0
    #     error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, i-1, -1.0);
    #     @assert error == 0
    #     error = GRBsetdblattrelement(model, GRB_DBL_ATTR_UB, i-1, 1.0);
    #     @assert error == 0
    #   end
    # elseif oa_iter == 3
    #   for i in 1:n
    #     @assert error == 0
    #     # error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, i-1, (vx[i] < 0 ? vx[i] * 1.5 : vx[i] / 1.5));
    #     error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, i-1, PI.xlo[i]);
    #     @assert error == 0
    #     error = GRBsetdblattrelement(model, GRB_DBL_ATTR_UB, i-1, PI.xhi[i]);
    #     @assert error == 0
    #   end
    # end
    
    # update active-z bounds and constraints
    for l in 1:L
      for i in 1:m[l]
        idx = n+sum(m[1:l-1])+Cint(i)
        if active[l][i] && !active_prev[l][i]
          # zbounds
          # println("----------\nl=$l, i=$i, idx=$idx, (x,z,t)idx = $((xidx,zidx,tidx))\n----------")
          F64P_lb = Ref{Float64}(NaN)
          F64P_ub = Ref{Float64}(NaN)
          error = GRBgetdblattrelement(model, GRB_DBL_ATTR_LB, idx-1, F64P_lb);
          @assert error == 0
          error = GRBgetdblattrelement(model, GRB_DBL_ATTR_UB, idx-1, F64P_ub);
          @assert error == 0
          error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, idx-1, 0.0);
          @assert error == 0
          error = GRBsetdblattrelement(model, GRB_DBL_ATTR_UB, idx-1, Cdouble(+Inf));
          @assert error == 0
          # println("updating bounds from ($(F64P_lb[]), $(F64P_ub[])) to (0,+Inf) for z[l][i]")
          
          # zconstraint
          #  (Ax)_W - z_W - t1_W <= -bW + (r/k)1_W
          # println("adding constraint A[l][i,:]' * x - z[l][i] - t[l] <= -(b[l][i] - r[l]/k[l])") # per BasRockRoy(eq20)
          if issparse(PI.A[l])
            cval = [PI.A[l][i,:].nzval; -1.0; -1.0]
            cind = [PI.A[l][i,:].nzind; idx; n+M+l]
            cnnz = length(cval)
          else
            cval = [PI.A[l][i,:]; -1.0; -1.0]
            cind = [collect(1:n); idx; n+M+l]
            cnnz = length(cval)
          end
          error = GRBaddconstr(
            model,              # GRBmodel   : *model
            Base.Csize_t(cnnz), # size_t     : numnz
            Cint.(cind.-1),     # int        : *cind
            Cdouble.(cval),     # double     : *cval
            GRB_LESS_EQUAL,     # char       : sense
            rhs_c2[l][i],       # double     : *rhs
            C_NULL,             # const char : **constrname
          )
          @assert error == 0
          order_added += 1
          condict[(i,l)] = order_added
          zdict[(i,l)] = idx-n
        else
          F64P_lb = Ref{Float64}(NaN)
          F64P_ub = Ref{Float64}(NaN)
          error = GRBgetdblattrelement(model, GRB_DBL_ATTR_LB, idx, F64P_lb);
          @assert error == 0
          error = GRBgetdblattrelement(model, GRB_DBL_ATTR_UB, idx, F64P_ub);
          @assert error == 0
          # println("z[l=$l][i=$i] bounds are ($(F64P_lb[]), $(F64P_ub[]))")
        end
      end
    end
    error = GRBupdatemodel(model)
    @assert error == 0

    # solve subproblem
    try
      solve_time = @elapsed error = GRBoptimize(model)
      error = GRBgetdblattr(model, "RunTime", RunTime);
      @assert error == 0
      res[:times][oa_iter] = RunTime[]
      if method == 2
        error = GRBgetintattr(model, "BarIterCount", BarIterCount);
        @assert error == 0
        res[:iters][oa_iter] = BarIterCount[]
      else
        error = GRBgetdblattr(model, "IterCount", IterCount);
        @assert error == 0
        res[:iters][oa_iter] = IterCount[]
      end
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
    
    # check stop
    stop = true
    for l in 1:L
      stop *= all(active[l] .== active_prev[l])
    end
    if stop || oa_iter > options[:max_oa_iter]
      if tolsub == min(options[:relgap_tol_hi], options[:pinfeas_tol_hi], options[:dinfeas_tol_hi])
        retcode = -1
        break
      else
        # solve with desired tolerance
        tolsub = min(options[:relgap_tol_hi], options[:pinfeas_tol_hi], options[:dinfeas_tol_hi])
        error = GRBupdatemodel(model)
        @assert error == 0
        error = GRBsetdblparam(GRBgetenv(model), "BarConvTol", tolsub);
        @assert error == 0
        error = GRBsetdblparam(GRBgetenv(model), "FeasibilityTol", tolsub);
        @assert error == 0
        error = GRBsetdblparam(GRBgetenv(model), "OptimalityTol", tolsub);
        @assert error == 0
        error = GRBupdatemodel(model)
        @assert error == 0
        solve_time = @elapsed error = GRBoptimize(model)
        error = GRBgetdblattr(model, "RunTime", RunTime);
        @assert error == 0
        res[:times][oa_iter+1] = RunTime[]
        if method == 2
          error = GRBgetintattr(model, "BarIterCount", BarIterCount);
          @assert error == 0
          res[:iters][oa_iter+1] = BarIterCount[]
        else
          error = GRBgetdblattr(model, "IterCount", IterCount);
          @assert error == 0
          res[:iters][oa_iter+1] = IterCount[]
        end
      end
    elseif time() - start_time > options[:maxtime]
      retcode = -2
      break
    else
      for l in 1:L
        active_prev[l] .= active[l]
      end
    end
  end
  
  #
  # recover solution
  #

  NumVars = Ref{Cint}();
  NumConstrs = Ref{Cint}();
  ObjVal = Ref{Cdouble}(); # primal obj
  ObjBound = Ref{Cdouble}(); # dual obj
  xzt = zeros(Cdouble, nn); # xzt
  # nu_2 = [zeros(Cdouble, sum(active[l])) for l in 1:L]; # c2
  nu_2 = zeros(Cdouble, sum(sum(active)))
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
  error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, nn, xzt);
  @assert error == 0
  error = GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, L, sum(sum(active)), nu_2);
  # idx_lo = idx_hi = 0
  # for l in 1:L
  #   if l == 1
  #     idx_lo = 0 + L
  #   else
  #     idx_lo = sum(sum(active[k]) for k in 1:l-1) + L
  #   end
  #   len = sum(sum(active[l]))
  #   idx_hi = idx_lo + len-1 + L
  #   error = GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, idx_lo, len, nu_2[l]);
  #   @assert error == 0
  # end
  error = GRBgetdblattrarray(model, GRB_DBL_ATTR_RC, 0, n, mu);
  @assert error == 0
  error = GRBgetdblattr(model, "DNumNZs", DNumNZs);
  @assert error == 0

  # AA = zeros(mm, nn)
  # for i in 1:mm
  #   for j in 1:nn
  #     vv = Ref{Float64}(NaN)
  #     error = GRBgetcoeff(model, Cint(i-1), Cint(j-1), vv)
  #     @assert error == 0
  #     AA[i,j] = vv[]
  #   end
  # end
  # bdbd_lo = zeros(Cdouble, nn);
  # bdbd_hi = zeros(Cdouble, nn);
  # error = GRBgetdblattrarray(model, GRB_DBL_ATTR_LB, 0, nn, bdbd_lo);
  # @assert error == 0
  # error = GRBgetdblattrarray(model, GRB_DBL_ATTR_UB, 0, nn, bdbd_hi);
  # @assert error == 0

  
  # populate result
  # res[:conmat] = AA
  # res[:bds] = hcat(bdbd_lo, bdbd_hi)
  res[:nnz] = DNumNZs[]
  res[:numvar] = NumVars[]
  res[:numcon] = NumConstrs[]
  res[:walltime] = time() - start_time
  res[:pobj] = ObjVal[]
  res[:x] = xzt[xidx]
  res[:normx] = norm(res[:x])
  res[:zz] = xzt[zidx]
  res[:z] = [zeros(m[l]) for l in 1:L]
  res[:t] = xzt[tidx]
  res[:mu] = mu
  res[:nu] = nu_2
  res[:condict] = condict
  res[:zdict] = zdict
  res[:active] = active
  res[:lambda] = [zeros(m[l]) for l in 1:L]
  for l in 1:L
    for i in 1:PI.m[l]
      if haskey(res[:condict], (i,l))
        conidx = res[:condict][(i,l)]
        res[:lambda][l][i] = res[:nu][conidx]
      end
      if haskey(res[:zdict], (i,l))
        zidx = res[:zdict][(i,l)]
        res[:z][l][i] = res[:zz][zidx]
      end
    end
  end
  res[:lambda] = vcat(res[:lambda]...)
  res[:z] = vcat(res[:z]...)
  res[:final_subproblem_status] = (Status[] == 2 ? 1 : 0)
  res[:status] = retcode > 0 ? 1 : 0
  res[:retcode] = retcode
  res[:times] = res[:times][1:oa_iter]
  res[:iters] = res[:iters][1:oa_iter]
  res[:iter] = length(findall(res[:iters] .> 0))
  # res[:y] = vcat([PI.A[l] * res[:x] .+ PI.b[l] for l in 1:L]...)
  
  #
  # free model
  #
  
  @assert error == 0
  error = GRBfreemodel(model)
  @assert error == 0
  GRBfreeenv(env)

  # clear memory
  Acsr = nothing; GC.gc()
  println("collected Acsr constraint")
  
  return res
end