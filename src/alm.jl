function alm(
  x0::Vf, lambda0::Vector{Vf}, mux0::Vf, sigma0::Tf, tau0::Tf, PI::ProblemInstance,
  AP0::ALMParameters, SP0::SSNParameters, kgrid::Union{Nothing,AbstractVector}, debug::Bool=false, verb::Int64=1, fullit::Bool=false, zero_x::Bool=false
) where {Tf<:AbstractFloat, Vf<:AbstractVector{Tf}}
  alm_start = time();
  SP = deepcopy(SP0)
  AP = deepcopy(AP0)
  # handle gridsolve
  if isnothing(kgrid)
    kgrid = [PI.k]
  end
  res_grid = Dict()
  # start grid clock
  grid_start = time();
  warm = false
  BXratio = AP.BXratio
  println("ALM: BXratio=$BXratio")

  #
  # initialize
  #
  
  Ti = Int64;
  m, n, L = PI.m, PI.n, PI.L;
  
  # primal: x
  x = zeros(Tf, n);
  xssn = zeros(Tf, n);
  xls = zeros(Tf, n);
  xhat = zeros(Tf, n);
  xprev = zeros(Tf, n);
  xbest = zeros(Tf, n);
  xtest = zeros(Tf, n);
  x .= x0;
  xssn .= x0;
  xprev .= x0;
  
  # primal: y
  ybar = [zeros(Tf, m[l]) for l in eachindex(m)];
  ybest = [zeros(Tf, m[l]) for l in eachindex(m)];
  ytest = [zeros(Tf, m[l]) for l in eachindex(m)];
  ytilde = [zeros(Tf, m[l]) for l in eachindex(m)];
  yhat = [zeros(Tf, m[l]) for l in eachindex(m)];
  ytmp = [zeros(Tf, PI.m[l]) for l in eachindex(PI.m)];
  y = ybar; #! pointer
  yprev = [zeros(Tf, m[l]) for l in eachindex(m)];
  
  # primal: y helpers
  sig = [zeros(Int64, m[l]) for l in eachindex(m)]; # permutation
  k0bar = zeros(Ti, PI.L);
  k1bar = zeros(Ti, PI.L);
  
  # primal: u
  ubar::Vf = zeros(Tf, PI.n);
  ubest::Vf = zeros(Tf, PI.n);
  uhat::Vf = zeros(Tf, PI.n);
  utilde::Vf = zeros(Tf, PI.n);
  u = ubar; #! pointer
  uprev::Vf = zeros(Tf, PI.n);
  
  # dual: lambda
  lambda = [zeros(Tf, m[l]) for l in eachindex(m)];
  lambdabest = [zeros(Tf, m[l]) for l in eachindex(m)];
  for l in 1:PI.L
    lambda[l] .= lambda0[l];
  end
  sinvlambda = [zeros(Tf, m[l]) for l in eachindex(m)];
  lambdaprev = [zeros(Tf, m[l]) for l in eachindex(m)];
  lambdahat = [zeros(Tf, m[l]) for l in eachindex(m)];
  
  # dual: mu
  mux::Vf = zeros(Tf, PI.n);
  muxbest::Vf = zeros(Tf, PI.n);
  mux .= mux0;
  muxprev::Vf = zeros(Tf, PI.n);
  muxhat::Vf = zeros(Tf, PI.n);
  sinvmux::Vf = zeros(Tf, PI.n);
  negmux_neg = zeros(Tf, PI.n); # [-mux]^-
  negmux_pos = zeros(Tf, PI.n); # [-mux]^+
  
  # helper: newton direction
  d = zeros(Tf, n);
  diagDinv::Vf = zeros(Tf, n);
  diagDinvsqrt::Vf = zeros(Tf, n);
  diagDsqrt::Vf = zeros(Tf, n);
  Dinv_g::Vf = zeros(Tf, n);
  Ttilde = zeros(Tf, sum(PI.m)+PI.L, PI.n);
  Ttilde_rowmax = Base.RefValue{Ti}(0);
  sol_M_swm = zeros(Tf, sum(PI.m)+PI.L);
  sol_n_swm = zeros(Tf, PI.n);
  c = [zeros(Rational{Int64}, 3) for l in eachindex(m)];
  a = [zeros(Tf, PI.n) for l in eachindex(m)];
  b = [zeros(Tf, PI.n) for l in eachindex(m)];
  active = zeros(Bool, PI.L); # active CVaR constraints
  mks = zeros(Tf, PI.L);
  diagjX = zeros(Tf, PI.n); # diagonal of jacobian of projection onto X
  
  # helper: matvec
  Ax = [zeros(Tf, m[l]) for l in eachindex(m)];
  Ad = [zeros(Tf, m[l]) for l in eachindex(m)];
  Atlampmu = zeros(Tf, n); # A' * lambda + mu
  Btlam = [zeros(Tf, m[l]) for l in eachindex(m)]; # Atildeinv' * lambda
  Cx = zeros(Tf, n);

  # objective
  obj::FunctionValue = FunctionValue(Tf, PI.n, PI.objtype) # 2==diagonal
  if !isnothing(PI.C)
    obj.H .= PI.C
  else
    obj.g .= PI.c
  end

  # max memory
  if PI.n >= AP.maxn || PI.n >= 10sum(PI.m)*sum(PI.k)/sum(PI.n)
    phi_full = FunctionValue(Tf, n, PI.objtype) # phi.H not used
  else
    phi_full = FunctionValue(Tf, PI.n, 3) # phi.H used
  end
  phi_smw = FunctionValue(Tf, n, PI.objtype) # phi.H not used
  phi = Base.Ref{FunctionValue}(phi_smw) #! pointer
  
  # diagnostics
  Ry = [zeros(Tf, m[l]) for l in eachindex(m)];
  Ru = zeros(Tf, PI.n);
  Rd = zeros(Tf, PI.n);

  # parameters
  sigma::Tf = one(Tf);
  sigma = sigma0;
  sigma_prev = zeros(Tf, 7)
  tau::Tf = zero(Tf);
  tau = tau0;
  stoptime = 7200
  maxnit0 = SP.maxnit
  if PI.objtype == 2
    if sigma0 < 0
      if PI.n <= 100
        SP.maxnit = 1000
        sigma = 5e-2
        tau = 0.0
      elseif PI.n <= 1_000
        SP.maxnit = 1000
        sigma = 5e-2
        tau = 0.0
      elseif PI.n > 1_000 && PI.n <= 10_000
        SP.maxnit = 1000
        sigma = 5e-1
        tau = 0.0
      elseif PI.n > 10_000 && PI.n <= 100_000
        SP.maxnit = 1000
        sigma = 1.0
        tau = 0.0
      elseif PI.n > 100_000
        SP.maxnit = 1000
        sigma = 5.0
        tau = 0.0
      end
    else
      sigma = sigma0
    end
  elseif PI.objtype == 1
    if sigma0 < 0 && tau < 0
      SP.maxnit = 100
      sigma = 1.0
      tau = 10.0
    elseif sigma0 < 0 && tau0 > 0
      SP.maxnit = 100
      sigma = 1.0
      tau = tau0
    elseif sigma0 > 0 && tau0 < 0
      SP.maxnit = 100
      sigma = sigma0
      tau = 1.0
    end
  end
  if verb >= 3
    println("(sigma,tau) = $((sigma,tau))")
  end

  #
  # start gridsolve
  #

  for (ii, kgrid_ii) in enumerate(kgrid)
    # start alm clock
    if length(kgrid) > 1
      alm_start = time(); # amortize initial allocation cost
    end
    PI.k .= kgrid_ii
    PI.Bt .= Function[Atmatinv_action(PI.m[l], PI.k[l]) for l in 1:PI.L]
    xbest .= xbest
    ubest .= ubest
    mux .= muxbest
    for l in 1:PI.L
      y[l] .= ybest[l]
      ybar[l] .= ybest[l]
      lambda[l] .= lambdabest[l]
    end
    println("\n\n\n     $ii, $kgrid_ii, $(PI.k)     \n\n\n")
    if ii >= 1
      warm = true
      if zero_x
        x .= 0.0
      end
    end
    
    # sigma/tau updates
    if warm
      sigma = (sigma0 * sigma)^(1/2)
      tau = (tau0 * tau)^(1/2)
    else
      sigma = sigma0
      tau = tau0
    end
    sigmamax = 1e6;
    sigmamin = 1e-9;
    taumax = 1e3;
    taumin = 1e-9;
    pinfeas_target_a = 1e-2; # phase 1+ pinfeas goal
    pinfeas_target = pinfeas_target_a;
    pinfeas_target_b = 1e-5;
    dinfeas_target = 1e-4; # phase 2+ pinfeas goal
    dinfeas_target_b = 1e-5
    phase = 1
    last_phase = 1
    iter_target = 15
    SP.maxnit = iter_target
    iter_target_satisfied = true
    phase_switch_min = 3
    phase_counter = 0
    good_subproblem_term = 0
    bad_subproblem_term = 0
    last_sigma_update = 0
    last_tau_update = 0
    aggressive = false

    # sort buffer
    sort_buffer = 0.2 # partial sort of min{m, |α_prev+β_prev| * (1 + sort_buffer) } entries for new |α+β|
    redo = [true for l in 1:PI.L] # re-do the partial sort?
    
    # tracking
    it = zeros(Tf, 3+9+6+3+1+1+1, AP.maxait);
    #=
    for each alm iteration, track:
      (3) -> ssn_iter, ssn_lsiter, ssn_time
      (9) -> ssn_time_nd, ssn_time_H, ssn_time_g, ssn_time_v, ssn_time_Ax, ssn_err, ssn_status, time_proj, time_sort
      (6) -> num_sps, time_Atlam, time_diagnostics, avg_nd_size, avg_step_size, partial_sort_count
      (3) -> ssn_pinfeas, ssn_dinfeas, ssn_relgap
      (1) -> time_assign
      (1) -> absbeta
      (1) -> absalpha
    =#

    sps::Ti = zero(Ti); # total strict CVaR projections
    pwin::Ti = 0; # primal wins
    dwin::Ti = 0; # dual wins
    err::Tf = +Inf;
    errbest::Tf = +Inf;
    # maxnit0 = SP.maxnit;
    tolsub0 = SP.tolsub;
    res = Dict();
    res[:walltime_i1] = NaN
    res[:walltime_i2] = NaN
    res[:walltime] = NaN
    almit_stop = 0.0;
    ctime_Atlam = 0.0
    subproblem_info = Dict(
      :start_time => alm_start,
      :time => 0.0, # time solving subproblem
      :iter => 0, # subproblem iterations
      :sp => 0, # strict projections
      :lsiter => 0, # subproblem ls iterations
      :ctime_proj => 0.0, # cumulative projection time
      :ctime_assign => 0.0, # cumulative assignment time
      :ctime_sort => 0.0, # cumulative projection time
      :ctime_nd => 0.0, # cumulative newton solve time
      :ctime_Hphi => 0.0, # cumulative jac time
      :ctime_gphi => 0.0, # cumulative grad time
      :ctime_vphi => 0.0, # cumulative obj time
      :ctime_Ax => 0.0, # cumulative A * x time (via A * d)
      :ctime_diagnostics => 0.0, # cumulative diagnostic time
      :err => Inf, # ||∇φ||
      :status => 0,
      :avg_nd_size => 0, # |β+1|
      :absbeta => 0, # |β+1|
      :absalpha => 0, # |β+1|
      :avg_step_size => 0.0,
      :partial_sort_count => 0,
      :normg => 0.0,
      :normg_prev => 0.0,
      :pinfeas => NaN,
      :dinfeas => NaN,
      :relgap => NaN,
      :pobj => NaN,
      :dobj => NaN,
    )

    # iterations
    ait::Ti = zero(Ti);
    almstop::Bool = false;
    retcodealm::Ti = 0;

    # diagnostics
    pinfeas, y_Axpb, u_x, dinfeas, relgap, objp, objd = +Inf, +Inf, +Inf, +Inf, +Inf, +Inf, +Inf;

    #
    # iteration 0
    #
    
    # get y(x) and Ax
    println("getting sorted/projected y(x) = Proj(Ax + b - sinvlambda)")
    subproblem_info[:ctime_Ax] += @elapsed begin
    for l in 1:PI.L # Threads.@threads 
      sinvlambda[l] .= lambda[l] ./ sigma
      # compute via BLAS: Ax[l] .= PI.A[l] * x
      # BLAS: gemv!(tA, alpha, A, x, beta, y) => y = alpha*A*x + beta*y
      # @views BLAS.gemv!('N', true, PI.A[l], x, false, Ax[l])
      @views mul!(Ax[l], PI.A[l], x, true, false)
      ytilde[l] .= Ax[l] .+ PI.b[l] .- sinvlambda[l]
    end; end # ctime_Ax
    subproblem_info[:ctime_sort] += @elapsed begin
    for l in 1:PI.L
      if PI.k[l] >= 1
        # @views sortperm!(sig[l], ytilde[l], rev=true)
        @views sortperm!(sig[l], ytilde[l], rev=true, alg=PartialQuickSort(1:PI.k[l]))
      else
        @views sortperm!(sig[l], ytilde[l], rev=true, alg=PartialQuickSort(1:PI.k[l]))
      end
      @views mks[l] = sum(ytilde[l][sig[l][1:PI.k[l]]])
      active[l] = mks[l] > PI.r[l]
    end; end # ctime_sort
    subproblem_info[:ctime_proj] += @elapsed begin
    redo .= true
    sp = project_B!(ytilde, sig, ybar, yhat, ytmp, k0bar, k1bar, PI.r, PI.k, active, redo)
    end # ctime_proj
    sps += sp # separated out of subproblem
    
    # get u(x)
    subproblem_info[:ctime_assign] += @elapsed begin
    sinvmux .= mux ./ (BXratio * sigma);
    utilde .= x .- sinvmux;
    end # ctime_assign
    subproblem_info[:ctime_proj] += @elapsed begin
    project_X!(utilde, ubar, uhat, PI.xlo, PI.xhi)
    end # ctime_proj

    # get phi
    subproblem_info[:ctime_vphi] += @elapsed begin
    # vgphi!(phi[], obj, x, x0, yhat, uhat, sigma, tau, PI, Cx, BXratio)
    vgphi!(phi[], obj, x, xprev, yhat, uhat, sigma, tau, PI, Cx, BXratio)
    #! phi_smw is in use at this point
    phi_full.v = phi[].v; phi_smw.v = phi[].v
    phi_full.g .= phi[].g; phi_smw.g .= phi[].g
    phi_full.d = phi[].d; phi_smw.d = phi[].d
    
    end # ctime_vphi

    # get Atlam + mu, Btlam
    ctime_Atlam += @elapsed begin
    Atlampmu .= mux
    for l in 1:PI.L # Threads.@threads
      @views mul!(Atlampmu, PI.A[l]', lambda[l], true, true)
    end; end # ctime_Atlam
    
    # diagnostics
    subproblem_info[:ctime_diagnostics] += @elapsed begin
    (pinfeas, y_Axpb, u_x, dinfeas, relgap, objp, objd, _) = calc_ALM_line(
      x, xprev, ybar, ubar, lambda, mux, negmux_neg, negmux_pos, Ax, Btlam, Atlampmu, phi[], obj, sigma, tau, PI,
      Cx, Ry, Ru, Rd, sig
    )
    end # ctime_diagnostics
    err = max(abs(pinfeas), abs(dinfeas), abs(relgap))
    if err <= errbest
      subproblem_info[:ctime_assign] += @elapsed begin
      xbest .= x;
      ybest .= y;
      ubest .= u;
      lambdabest .= lambda;
      muxbest .= mux;
      errbest = err;
      end # ctime_assign
    end
      
    # print
    subproblem_info[:ctime_diagnostics] += @elapsed begin
    if verb >= 1
      print(ALM_line(ait+1/2, pinfeas, pinfeas_target, dinfeas, dinfeas_target, relgap, objp, objd, sigma, tau, subproblem_info[:lsiter], subproblem_info[:iter], phase, SP.tolsub, subproblem_info[:normg], 0.0, 0.0, true))
    end; end # ctime_diagnostics

    if false # err <= AP.stoptol
      println("subproblem found optimal solution: err=$err")
      subproblem_info[:ctime_assign] += @elapsed begin
      res[:retcode] = retcodealm
      res[:status] = (retcodealm != 4 ? 1 : 0) # 1=solved, 0=notsolved
      res[:x] = deepcopy(xbest)
      res[:y] = deepcopy(ybest)
      res[:u] = deepcopy(ubest)
      res[:lambda] = deepcopy(lambdabest)
      res[:mu] = deepcopy(muxbest)
      res[:pinfeas] = pinfeas
      res[:dinfeas] = dinfeas
      res[:relgap] = relgap
      res[:pobj] = objp
      res[:dobj] = objd
      res[:cvar] =  [sum( (PI.A[l] * xbest .+ PI.b[l])[sig[l][1:PI.k[l]]] ) / PI.k[l] for l in 1:PI.L]
      # res[:maxksum] = copy([maxksum(ybar[l], PI.k[l]) for l in 1:PI.L])
      res[:maxksum] = res[:cvar] .* PI.k
      res[:maxksumrhs] = PI.r
      res[:k] = kgrid_ii
      res[:walltime] = time()-alm_start
      res[:tau] = tau
      res[:sigma] = sigma
      res[:iter] = 0
      res_grid[ii] = res
      end # ctime_assign
    else

      #
      # perform alm iters
      #
      
      while true

        ait += 1
        almit_start = time()
        
        #
        # ait: subproblem tolerance
        #

        subproblem_info[:ctime_assign] += @elapsed begin
        xhat .= x .- xprev;
        yhat .= y .- yprev;
        uhat .= u .- uprev;
        lambdahat .= lambda .- lambdaprev;
        muxhat .= mux .- muxprev;
        end # ctime_assign
        subproblem_info[:ctime_assign] += @elapsed begin
        if tau > 0.0
          const1 = min(sqrt(tau), 10.0) / sigma
          epsilon = 10.0 / (ait^(1.01))
          delta = min(epsilon, 0.95)
          tolsub_A = const1 * epsilon
          if ait <= 1
            tolsub_B = 1.0 # first subproblem
          else
            pdnorm = sqrt(tau * norm(xhat, 2)^2 + norm(yhat, 2)^2 + norm(uhat, 2)^2 + norm(lambdahat, 2)^2 + norm(muxhat, 2)^2)
            tolsub_B = const1 * delta * pdnorm
          end

          # control number of ssn iterations and dinfeas
          if !warm && (PI.m > 50*PI.n)
            if (ait <= 5 && !warm) || (pinfeas <= 1e-12 && dinfeas >= 1e-3)
              SP.tolsub = 1.0 / ait^1.01
            elseif pinfeas > 1e-2 # (Inf, 1e-2)
              if dinfeas <= pinfeas/1.5 # good dinf
                SP.tolsub = max(min(10.0 / (ait^1.01), pinfeas), AP.stoptol/10)
              else
                SP.tolsub = max(min(10.0 / (ait^1.01), pinfeas/1.1), AP.stoptol/10)
              end
            elseif pinfeas > 1000AP.stoptol # (1e-2, 1e-3)
              if dinfeas <= pinfeas/2 # good dinf
                SP.tolsub = max(min(1.5dinfeas, pinfeas, tolsub_A, tolsub_B), AP.stoptol/10)
              else
                SP.tolsub = max(min(1.5dinfeas, pinfeas/1.2, tolsub_A, tolsub_B), AP.stoptol/10)
              end
            elseif pinfeas > 50AP.stoptol # (1e-3, 5e-5)
              if dinfeas <= pinfeas/2.5 # good dinf
                SP.tolsub = max(min(1.2dinfeas, pinfeas, tolsub_A, tolsub_B), AP.stoptol/1.1)
              else
                SP.tolsub = max(min(1.2dinfeas, pinfeas/1.5, tolsub_A, tolsub_B), AP.stoptol/1.1)
              end
            else # (5e-5, 0.0)
              if dinfeas <= pinfeas/2.5 # good dinf
                SP.tolsub = max(min(1.1dinfeas, pinfeas, tolsub_A, tolsub_B), AP.stoptol/1.5)
              else
                SP.tolsub = max(min(1.1dinfeas, pinfeas/1.6, tolsub_A, tolsub_B), AP.stoptol/1.5)
              end
            end
          else
            # println("~here~")
            # warmstart
            SP.tolsub = max(min(1e-6, pinfeas, dinfeas), (AP.stoptol <= 1e-8 ? AP.stoptol/1.1 : AP.stoptol/5))
          end
        else
          SP.tolsub = max(min(pinfeas / 2.0, 1e-1), AP.stoptol)
          if verb >= 3
            println(@sprintf("tolsub = %0.2e", SP.tolsub))
          end
        end
        
        # if PI.n > sum(PI.m) || fullit
        if 2PI.n >= sum(PI.m)/sum(PI.L) || fullit          
          SP.maxnit = maxnit0
          println("SP.maxnit  =  $(SP.maxnit)")
        else
          if pinfeas >= dinfeas/10000 || ait <= 5 # if pinfeas is not too much better than dinfeas
            if (pinfeas > pinfeas_target_a) || (ait < 5 && !warm) # if pinfeas in [Inf, 1e-3], afford fewer iterations
              if subproblem_info[:normg] <= subproblem_info[:normg_prev] && phase == 1
                SP.maxnit = min(2iter_target, max(iter_target, SP.maxnit-5)) # between 15 and 30
                # println("ssnit = $(SP.maxnit): (1)")
              elseif subproblem_info[:normg] <= dinfeas
                SP.maxnit = min(2iter_target, max(iter_target, SP.maxnit+3)) # between 15 and 30
              else
                SP.maxnit = min(4iter_target, max(iter_target, SP.maxnit+3)) # between 15 and 60
                # println("ssnit = $(SP.maxnit): (2)")
              end
            elseif (pinfeas > pinfeas_target_b) # if pinfeas in [1e-3, 1e-5]
              if subproblem_info[:normg] <= subproblem_info[:normg_prev]/2 && phase == 1
                SP.maxnit = min(4iter_target, max(iter_target, SP.maxnit+5)) # between 15 and 60
                # println("ssnit = $(SP.maxnit): (3)")
              else
                SP.maxnit = min(4iter_target, SP.maxnit+10)
                # println("ssnit = $(SP.maxnit): (4)")
              end
            else
              SP.maxnit = maxnit0
              # println("ssnit = $(SP.maxnit): (5)")
            end
          else
            SP.maxnit = maxnit0
            # println("ssnit = $(SP.maxnit): (6)")
          end
        end
        println("ssnit = $(SP.maxnit)")
        end # ctime_assign

        #
        # ait: solve subproblem (update primal variables)
        #
        
        subproblem_info[:ctime_assign] += @elapsed begin
        xprev .= x;
        uprev .= u;
        yprev .= y;
        subproblem_info[:normg_prev] = subproblem_info[:normg];
        end # ctime_assign
        # ssn calling sequence: ssn!(x0 ≅ xlm, x ≅ xssn, ...)
        try
          if subproblem_info[:status] != 6
            # noncycling ssn
            ssn!(
              subproblem_info,
              x, xssn, xls, xhat, d, Ax, Ad,
              ytilde, ybar, sig,
              utilde, ubar,
              yhat, ytmp, k0bar, k1bar,
              uhat,
              diagDinv, diagDinvsqrt, diagDsqrt, Dinv_g, sol_M_swm, sol_n_swm,
              lambda, mux, negmux_neg, negmux_pos,
              sinvlambda, sinvmux,
              Atlampmu, Btlam,
              sigma, tau,
              phi_full, phi_smw, phi, obj,
              (pinfeas, y_Axpb, u_x, dinfeas, relgap, objp, objd),
              PI,
              Cx, Ry, Ru, Rd,
              Ttilde, Ttilde_rowmax, c, a, b, active, diagjX, mks,
              SP.maxnit, SP.maxlit, SP.tolsub, SP.ls_method, SP.ls_tau, SP.ls_c1, SP.ls_c2, SP.stagit, sort_buffer, SP.lin_solver,
              AP.stoptol, AP.maxn,
              debug, verb, BXratio
            )
          else
            # cycling ssn
            ssn!(
              subproblem_info,
              x, xssn, xls, xhat, d, Ax, Ad,
              ytilde, ybar, sig,
              utilde, ubar,
              yhat, ytmp, k0bar, k1bar,
              uhat,
              diagDinv, diagDinvsqrt, diagDsqrt, Dinv_g, sol_M_swm, sol_n_swm,
              lambda, mux, negmux_neg, negmux_pos,
              sinvlambda, sinvmux,
              Atlampmu, Btlam,
              sigma, tau,
              phi_full, phi_smw, phi, obj,
              (pinfeas, y_Axpb, u_x, dinfeas, relgap, objp, objd),
              PI,
              Cx, Ry, Ru, Rd,
              Ttilde, Ttilde_rowmax, c, a, b, active, diagjX, mks,
              SP.maxnit, SP.maxlit, SP.tolsub, SP.ls_method, cycle_ls_tau, cycle_ls_c1, cycle_ls_c2, SP.stagit, sort_buffer, SP.lin_solver,
              AP.stoptol, AP.maxn,
              debug, verb, BXratio
            )
          end
        catch e
          println(e)
          showerror(stdout, e, catch_backtrace())
          almstop = true
          retcodealm = -1
        end
        subproblem_info[:ctime_assign] += @elapsed begin
        x .= xssn # update xalm
        sps += subproblem_info[:sp] # strict projections
        if subproblem_info[:status] == 1 || subproblem_info[:status] == 2
          good_subproblem_term += 1
          bad_subproblem_term = 0
        else
          good_subproblem_term = 0
          bad_subproblem_term += 1
        end
        end # ctime_assign
        
        #
        # ait: subproblem stopping
        #

        if subproblem_info[:err] <= AP.stoptol
          # solved in subproblem
          almstop = true
          retcodealm = 1
          pinfeas, dinfeas, relgap = subproblem_info[:pinfeas], subproblem_info[:dinfeas], subproblem_info[:relgap]
          pobj, dobj = subproblem_info[:pobj], subproblem_info[:dobj]
        end
        
        #
        # ait: dual (lambda, mu) updates
        #

        # update lambda
        subproblem_info[:ctime_assign] += @elapsed begin
        lambdaprev .= lambda
        for l in 1:L
          lambda[l] .+= sigma .* (ybar[l] .- (Ax[l] .+ PI.b[l]))
          lambda[l][abs.(lambda[l]) .< 1e-10] .= 0.0
        end
        
        # update mux
        muxprev .= mux
        mux .+= (sigma * BXratio) .* (ubar .- x)
        end # ctime_assign

        # get Atlam + mu, Btlam
        ctime_Atlam += @elapsed begin
        Atlampmu .= mux
        for l in 1:PI.L
          @views mul!(Atlampmu, PI.A[l]', lambda[l], true, true)
        end; end # ctime_Atlam
        
        #
        # ait: diagnostics
        #
        
        if retcodealm != 1
          subproblem_info[:ctime_diagnostics] += @elapsed begin
          (pinfeas, y_Axpb, u_x, dinfeas, relgap, objp, objd, _) = calc_ALM_line(
            x, xprev, y, u, lambda, mux, negmux_neg, negmux_pos, Ax, Btlam, Atlampmu, phi[], obj, sigma, tau, PI,
            Cx, Ry, Ru, Rd, sig
          )
          end # ctime_diagnostics
        end
        err = max(abs(pinfeas), abs(dinfeas), abs(relgap))
        if err <= errbest
          subproblem_info[:ctime_assign] += @elapsed begin
          xbest .= x
          ybest .= y
          ubest .= u
          lambdabest .= lambda
          muxbest .= mux
          errbest = err
          end # ctime_assign
        end
        
        # print
        almit_stop = time()
        subproblem_info[:ctime_diagnostics] += @elapsed begin
        if verb >= 1
          print(ALM_line(ait, pinfeas, pinfeas_target, dinfeas, dinfeas_target, relgap, objp, objd, sigma, tau, subproblem_info[:lsiter], subproblem_info[:iter], phase, SP.tolsub, subproblem_info[:normg], almit_stop-almit_start, time()-alm_start, verb >= 2))
        end; end # ctime_diagnostics

        #
        # ait: alm stopping
        #

        if err <= AP.stoptol_i1 && isnan(res[:walltime_i1]) # intermediate time
          res[:walltime_i1] = time()-alm_start
        end
        if err <= AP.stoptol_i2 && isnan(res[:walltime_i2]) # intermediate time
          res[:walltime_i2] = time()-alm_start
        end
        if err <= AP.stoptol
          # solved in ait (with dual update)
          almstop = true
          retcodealm = 2
        end
        
        #
        # ait: tracking
        #

        subproblem_info[:ctime_assign] += @elapsed begin
        it[1,  ait] = subproblem_info[:iter]; subproblem_info[:iter] = 0
        it[2,  ait] = subproblem_info[:lsiter]; subproblem_info[:lsiter] = 0
        it[3,  ait] = subproblem_info[:time]; subproblem_info[:time] = 0
        it[4,  ait] = subproblem_info[:ctime_nd]; subproblem_info[:ctime_nd] = 0
        it[5,  ait] = subproblem_info[:ctime_Hphi]; subproblem_info[:ctime_Hphi] = 0
        it[6,  ait] = subproblem_info[:ctime_gphi]; subproblem_info[:ctime_gphi] = 0
        it[7,  ait] = subproblem_info[:ctime_vphi]; subproblem_info[:ctime_vphi] = 0
        it[8,  ait] = subproblem_info[:ctime_Ax]; subproblem_info[:ctime_Ax] = 0
        it[9,  ait] = subproblem_info[:err]; subproblem_info[:err] = 0
        it[10, ait] = subproblem_info[:status]; subproblem_info[:status] = 0
        it[11, ait] = subproblem_info[:ctime_proj]; subproblem_info[:ctime_proj] = 0
        it[12, ait] = subproblem_info[:ctime_sort]; subproblem_info[:ctime_sort] = 0
        it[13, ait] = sps; sps = 0
        it[14, ait] = ctime_Atlam; ctime_Atlam = 0
        it[15, ait] = subproblem_info[:ctime_diagnostics]; subproblem_info[:ctime_diagnostics] = 0
        it[16, ait] = subproblem_info[:avg_nd_size]; subproblem_info[:avg_nd_size] = 0
        it[17, ait] = subproblem_info[:avg_step_size]; subproblem_info[:avg_step_size] = 0
        it[18, ait] = subproblem_info[:partial_sort_count]; subproblem_info[:partial_sort_count] = 0
        it[19, ait] = subproblem_info[:pinfeas];
        it[20, ait] = subproblem_info[:dinfeas];
        it[21, ait] = subproblem_info[:relgap];
        # stricter ls params when cycling 
        cycle_ls_tau = 0.6
        cycle_ls_c1 = 1e-3
        cycle_ls_c2 = 0.2
        end # ctime_assign
        it[22, ait] = subproblem_info[:ctime_assign]; subproblem_info[:ctime_assign] = 0.0
        it[23, ait] = subproblem_info[:absbeta]; subproblem_info[:absbeta] = 0
        it[24, ait] = subproblem_info[:absalpha]; subproblem_info[:absalpha] = 0
        # push!(active_idx, findall(abs.(y[1] .- PI.A[1] * x .+ PI.b[1]) .< 1e-6))

        #
        # break
        #

        # GC.gc()

        if almstop
          break
        end

        #
        # ait: update ALM sigma and tau
        #

        if pinfeas < dinfeas
          pwin = pwin+1
          dwin = 0
        else
          pwin = 0
          dwin = dwin+1
        end
        iter_target_satisfied = it[1, ait] <= iter_target-1
        
        # phase selection
        # #=
        time_tausigma = @elapsed begin
        if phase == 1 # improve pinfeas
          phase_counter += 1
          if iter_target_satisfied
            if dwin > max(1, 1.5pwin) && dinfeas <= pinfeas / 2.0 # have dinfeas to spare
              sigmascale = 1.9 # attempt to improve pinfeas aggressively
              tauscale = 1.2 # help numerics and offset dwin > pwin
              pwin = dwin = 0
              # println("p1: good iter | spare dinfeas | suberr $(subproblem_info[:normg]) | suberrprev $(subproblem_info[:normg_prev]) | ctr = $phase_counter")
            else
              sigmascale = 1.5
              tauscale = 0.9
              # println("p1: good iter | normal | suberr $(subproblem_info[:normg]) | suberrprev $(subproblem_info[:normg_prev]) | ctr = $phase_counter")
            end
          else # "too many" subproblem iterations
            if dwin > max(1, 1.5pwin) && dinfeas <= pinfeas / 2.0 # have dinfeas to spare
              sigmascale = 1.3
              tauscale = 1.6 # help numerics and offset dwin > pwin
              pwin = dwin = 0
              # println("p1: bad iter | spare dinfeas | suberr $(subproblem_info[:normg]) | suberrprev $(subproblem_info[:normg_prev]) | ctr = $phase_counter")
            else
              sigmascale = 1.3
              tauscale = 1.1 # help numerics and offset dwin > pwin
              # println("p1: bad iter | normal | suberr $(subproblem_info[:normg]) | suberrprev $(subproblem_info[:normg_prev]) | ctr = $phase_counter")
            end
          end
          if pinfeas <= pinfeas_target_b && dinfeas <= AP.stoptol/2 # clean up pinfeas
            phase = 1
            if subproblem_info[:iter] <= 2iter_target
              sigmascale = 1.9
              tauscale = 1.0
            else
              sigmascale = 1.6
              tauscale = 1.05
            end
          elseif pinfeas <= pinfeas_target_b && dinfeas <= dinfeas_target_b/2 # improve pinfeas aggressively
            phase = 1
            if subproblem_info[:iter] <= 2iter_target
              sigmascale = 1.6
              tauscale = 1.0
            else
              sigmascale = 1.4
              tauscale = 1.05
            end
          elseif (pinfeas <= pinfeas_target || pinfeas <= dinfeas / 20.0) && ait > 5 && phase_counter >= phase_switch_min
            phase = 2
            pinfeas_target = min(pinfeas_target, pinfeas)
            phase_counter = 0
            # println("SWITCH: 1=>2, (pinfeas,dinfeas) target = $((pinfeas_target,dinfeas_target))")
          end
        elseif phase == 2 # improve dinfeas and maintain pinfeas
          phase_counter += 1
          if pinfeas > pinfeas_target && dinfeas <= dinfeas_target && phase_counter >= phase_switch_min
            phase = 1
            tauscale = 1.0/1.1
            sigmascale = 1.1
            # println("SWITCH: 2=>1: p2: bad pinfeas | suberr $(subproblem_info[:normg]) | suberrprev $(subproblem_info[:normg_prev]) | ctr = $phase_counter")
          else
            if iter_target_satisfied
              tauscale = 1/1.9
              sigmascale = 1.0
              # println("p2: good iter | normal | suberr $(subproblem_info[:normg]) | suberrprev $(subproblem_info[:normg_prev]) | ctr = $phase_counter")
            else
              if pwin > max(1, 1.5dwin) && pinfeas < pinfeas_target/5 # have pinfeas to spare
                tauscale = 1/1.9
                sigmascale = 0.9
                # println("p2: bad iter | spare pinfeas | suberr $(subproblem_info[:normg]) | suberrprev $(subproblem_info[:normg_prev]) | ctr = $phase_counter")
                pwin = dwin = 0
              else
                tauscale = 1/1.5
                sigmascale = 1.0
                # println("p2: bad iter | normal | suberr $(subproblem_info[:normg]) | suberrprev $(subproblem_info[:normg_prev]) | ctr = $phase_counter")
              end
            end
          end
          if dinfeas <= dinfeas_target && phase_counter >= phase_switch_min
            phase = 1
            dinfeas_target = min(dinfeas_target, pinfeas_target/1.5)
            phase_counter = 0
            # println("SWITCH: 2=>1, (pinfeas,dinfeas) target = $((pinfeas_target,dinfeas_target))")
          end
        end
        end # time_tausigma
        subproblem_info[:ctime_assign] += time_tausigma
        # println("time tausigma = $(time_tausigma)")
        
        # bad subproblem override
        if subproblem_info[:status] == 4 || subproblem_info[:status] == 6
          println("override")
          tauscale = 1.5
          sigmascale = 1/1.5
        end

        # infeasible override
        if sigma == sigmamax && all(prod(sigma_prev .== sigmamax))
          retcodealm = -1
          almstop = true
        end

        # record update
        if PI.objtype == 1
          if min(pinfeas, dinfeas) <= min(AP.stoptol*200) && max(pinfeas, dinfeas) <= AP.stoptol*500 && sum(PI.m) <= 5e6
            if PI.n < maximum(PI.m)
              tau = max(tau/3, taumin)
              tauscale = 1.0
              if verb >= 2
                println("override: tau = tau/3")
              end
            elseif phase == 2
              tau = min(tau/2, 1e-4)
              tauscale = 1.0
              if verb >= 2
                println("override: tau = min(tau/2, 1e-4)")
              end
            end
          end
        end
        sigma_prev[1] = sigma
        sigma_prev[2:end] = sigma_prev[1:end-1]
        sigma = max(sigmamin, min(sigmamax, sigma * sigmascale));
        if PI.objtype == 1
          tau = max(taumin, min(taumax, tau * tauscale));
        else
          tau = max(taumin, min(taumax, tau * tauscale));
        end
        
        #
        # ait+1/2: primal (y,z) updates
        #

        if abs(dinfeas) <= AP.stoptol/1.1 && abs(relgap) <= AP.stoptol/1.1
          # println("TESTING: tol=$(AP.stoptol), dinfeas=$dinfeas, relgap=$relgap")
          sigma_tests = [10^6, sigma]
        else
          sigma_tests = [sigma]
        end
        for sigma_test in sigma_tests
          # get y(x)
          for l in 1:PI.L
            # Ax[l] .= PI.A[l] * x #! updated within linesearch
            subproblem_info[:ctime_assign] += @elapsed begin
            sinvlambda[l] .= lambda[l] ./ sigma_test
            ytilde[l] .= Ax[l] .+ PI.b[l] .- sinvlambda[l]
            end # ctime_assign
            subproblem_info[:ctime_sort] += @elapsed begin
            idx_guess = min(Ti(ceil(k1bar[l] * (1.0 + sort_buffer))), PI.m[l])
            if PI.k[l] >= 1
              @views sortperm!(sig[l], ytilde[l], alg=PartialQuickSort(1:idx_guess), rev=true);
            else
              @views sortperm!(sig[l], ytilde[l], rev=true, alg=PartialQuickSort(1:PI.k[l]))
            end
            end # ctime_sort
            @views mks[l] = sum(ytilde[l][sig[l][1:PI.k[l]]])
            active[l] = mks[l] > PI.r[l]
          end
          subproblem_info[:ctime_proj] += @elapsed begin
          redo .= true
          sp = project_B!(ytilde, sig, ybar, yhat, ytmp, k0bar, k1bar, PI.r, PI.k, active, redo)
          end # ctime_proj
          sps += sp

          # ensure correct projection from partial sort
          psc = 0 # partial sort counter
          while true
            redo .= false
            for l in 1:PI.L
              redo[l] = !isempty(findall(ytilde[l][sig[l][k1bar[l]+1:PI.m[l]]] .> ybar[l][sig[l][k1bar[l]]]));
            end
            if prod(redo) == false
              break
            end
            psc += 1
            for l in 1:PI.L
              idx_guess = min(Ti(ceil(k1bar[l] * (1.0 + sort_buffer) * 2^psc)), PI.m[l])
              if PI.k[l] >= 1
                subproblem_info[:ctime_sort] += @elapsed @views sortperm!(sig[l], ytilde[l], alg=PartialQuickSort(1:idx_guess), rev=true); # additional partial sort
              else
                subproblem_info[:ctime_sort] += @elapsed @views sortperm!(sig[l], ytilde[l], rev=true, alg=PartialQuickSort(1:PI.k[l]))
              end
            end
            subproblem_info[:ctime_proj] += @elapsed begin
            sps += project_B!(
              ytilde, sig, 
              ybar, yhat, ytmp,
              k0bar, k1bar, PI.r, PI.k, active, redo
            )
            end # time_proj
          end
          subproblem_info[:partial_sort_count] += psc

          # get u(x)
          sinvmux .= mux ./ (BXratio * sigma_test)
          utilde .= x .- sinvmux
          subproblem_info[:ctime_proj] += @elapsed begin
          project_X!(utilde, ubar, uhat, PI.xlo, PI.xhi)
          end # ctime_proj

          # update phi
          subproblem_info[:ctime_vphi], subproblem_info[:ctime_gphi] = vgphi!(phi[], obj, x, xprev, yhat, uhat, sigma, tau, PI, Cx, BXratio)
          phi_full.v = phi[].v; phi_smw.v = phi[].v
          phi_full.g .= phi[].g; phi_smw.g .= phi[].g
          phi_full.d = phi[].d; phi_smw.d = phi[].d
      
          #
          # ait + 1/2: diagnostics
          #
          
          @assert(all(y[1].==ybar[1]))
          subproblem_info[:ctime_diagnostics] += @elapsed begin
          (pinfeas, y_Axpb, u_x, dinfeas, relgap, objp, objd, _) = calc_ALM_line(
            x, xprev, y, u, lambda, mux, negmux_neg, negmux_pos, Ax, Btlam, Atlampmu, phi[], obj, sigma, tau, PI,
            Cx, Ry, Ru, Rd, sig
          )
          end # ctime_diagnostics
          err = max(abs(pinfeas), abs(dinfeas), abs(relgap))
          if err <= errbest
            subproblem_info[:ctime_assign] += @elapsed begin
            xbest .= x
            ybest .= y
            ubest .= u
            lambdabest .= lambda
            muxbest .= mux
            errbest = err
            end # ctime_assign
          end

          #
          # ait + 1/2: alm stopping
          #

          if err <= AP.stoptol
            # solved in ait + 1/2 (with dual update)
            almstop = true
            retcodealm = 3
            if verb >= 2
              println("stopping! sigma = $sigma")
            end
          end
          
          # print
          almit_stop = time()
          subproblem_info[:ctime_diagnostics] += @elapsed begin
          if verb >= 2
            print(ALM_line(ait+1/2, pinfeas, pinfeas_target, dinfeas, dinfeas_target, relgap, objp, objd, sigma, tau, subproblem_info[:lsiter], subproblem_info[:iter], phase, SP.tolsub, subproblem_info[:normg], almit_stop-almit_start, time()-alm_start, verb >= 2))
          end; end # ctime_diagnostics

          if almstop
            break
          end
        end

        #
        # ait + 1/2: alm stopping
        #

        if err <= AP.stoptol
          # solved in ait + 1/2 (with dual update)
          almstop = true
          retcodealm = 3
        end
        
        #
        # ait + 1: stopping
        #

        if ait + 1 > AP.maxait
          almstop = true
          retcodealm = 4
        elseif time() - alm_start >= AP.maxtime
          almstop = true
          retcodealm = 5
        end
        
        #
        # break
        #

        if almstop
          break
        end
      
      end # while true
      
      #
      # return
      #

      # stop alm clock
      alm_stop = time()

      # reset SP parameters
      SP.maxnit = maxnit0
      SP.tolsub = tolsub0
      
      res[:retcode] = retcodealm
      res[:status] = ((retcodealm == -1 || retcodealm == 4 || retcodealm == 5) ? 0 : 1) # 1=solved, 0=notsolved
      res[:x] = deepcopy(xbest)
      res[:normx] = norm(res[:x])
      res[:y] = deepcopy(ybest)
      res[:k0bar] = deepcopy(k0bar)
      res[:k1bar] = deepcopy(k1bar)
      res[:u] = deepcopy(ubest)
      res[:lambda] = deepcopy(lambdabest)
      res[:mu] = deepcopy(muxbest)
      res[:cvar] =  [sum( (PI.A[l] * xbest .+ PI.b[l])[sig[l][1:PI.k[l]]] ) / PI.k[l] for l in 1:PI.L]
      res[:cvar_rhs] = PI.r./PI.k
      res[:cvar_infeas] = maximum(max.(res[:cvar] .- res[:cvar_rhs], 0.0))
      res[:pinfeas] = max(pinfeas, res[:cvar_infeas])
      res[:dinfeas] = dinfeas
      res[:relgap] = relgap
      res[:pobj] = objp
      res[:dobj] = objd
      res[:maxksum] = res[:cvar] .* PI.k
      res[:maxksumrhs] = PI.r
      res[:tau] = copy(tau)
      res[:sigma] = copy(sigma)
      res[:n] = PI.n
      res[:M] = sum(PI.m)
      res[:walltime] = alm_stop - alm_start
      res[:iter] = ait
      res[:nit] = it[1, 1:ait]
      res[:lsit] = it[2, 1:ait]
      res[:timessn] = it[3, 1:ait]
      res[:ctime_nd] = it[4, 1:ait]
      res[:ctime_Hphi] = it[5, 1:ait]
      res[:ctime_gphi] = it[6, 1:ait]
      res[:ctime_vphi] = it[7, 1:ait]
      res[:ctime_Ax] = it[8, 1:ait]
      res[:err] = it[9, 1:ait]
      res[:substatus] = it[10, 1:ait]
      res[:ctime_proj] = it[11, 1:ait]
      res[:ctime_sort] = it[12, 1:ait]
      res[:sps] = it[13, 1:ait]
      res[:ctime_Atlam] = it[14, 1:ait]
      res[:ctime_diagnostics] = it[15, 1:ait]
      res[:avg_nd_size] = it[16, 1:ait]
      res[:avg_step_size] = it[17, ait]
      res[:partial_sort_count] = it[18, 1:ait]
      res[:sub_pinfeas] = it[19, 1:ait]
      res[:sub_dinfeas] = it[20, 1:ait]
      res[:sub_relgap] = it[21, 1:ait]
      res[:ctime_assign] = it[22, 1:ait]
      res[:absbeta] = it[23, 1:ait]
      res[:absalpha] = it[24, 1:ait]
      res[:pct_timed] = (
        sum(res[:ctime_nd]) + sum(res[:ctime_Hphi]) + sum(res[:ctime_gphi]) + sum(res[:ctime_vphi]) +
        sum(res[:ctime_assign]) + sum(res[:ctime_proj]) + sum(res[:ctime_sort]) + sum(res[:ctime_Atlam]) + sum(res[:ctime_Ax]) +
        sum(res[:ctime_diagnostics])
      ) / res[:walltime]
      res[:objseed] = PI.objseed
      res[:conseed] = PI.conseed
      res[:k] = kgrid_ii
      res_grid[ii] = deepcopy(res)
    end
  end
  res_grid[:walltime] = time()-grid_start
  if length(kgrid) == 1
    res_grid[1][:walltime] = res_grid[:walltime]
    return res_grid[1], sigma, tau
  else
    return res_grid
  end
end
