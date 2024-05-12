#=
  in-place
=#
function nd!(
  xssn::Vf, xalm::Vf, d::Vf,
  ytilde::Vector{Vf}, ybar::Vector{Vf}, ytmp::Vector{Vf}, yhat::Vector{Vf}, sig::Vector{Vi}, k0bar::Vi, k1bar::Vi,
  utilde::Vf, ubar::Vf, uhat::Vf,
  sigma::Tf, tau::Tf,
  phi::FV1, obj::FV2, PI::ProblemInstance,
  Cx::Union{Nothing,Vf},
  Ttilde::Mf, Ttilde_rowmax::Base.RefValue{Ti}, c::Vector{Vr}, a::Vector{Vf}, b::Vector{Vf}, active::Vector{Bool}, diagjX::Vf,
  debug::Bool=false, BXratio::Tf=1.0
) where {
  Tf<:AbstractFloat,Vf<:AbstractVector{Tf},Mf<:AbstractMatrix{Tf},
  Ti<:Integer,Vi<:AbstractVector{Ti},
  Tr<:Rational{Ti},Vr<:AbstractVector{Tr},
  FV1<:FunctionValue,FV2<:FunctionValue
}
  ctime_H = Hphi!(
    xssn, xalm, d,
    ytilde, ybar, ytmp, sig, k0bar, k1bar,
    utilde, ubar,
    sigma, tau,
    phi, obj, PI,
    Cx,
    Ttilde, Ttilde_rowmax, c, a, b, active, diagjX,
    debug, BXratio
  )
  ctime_nd = @elapsed d .= phi.H \ -phi.g
  return ctime_H, ctime_nd
end

function nd_swm!(
  xssn::Vf, xalm::Vf, d::Vf,
  ytilde::Vector{Vf}, ybar::Vector{Vf}, ytmp::Vector{Vf}, yhat::Vector{Vf}, sig::Vector{Vi}, k0bar::Vi, k1bar::Vi,
  utilde::Vf, ubar::Vf, uhat::Vf,
  diagDinv::Vf, diagDinvsqrt::Vf, diagDsqrt::Vf, Dinv_g::Vf, sol_M_swm::Vf, sol_n_swm::Vf, 
  sigma::Tf, tau::Tf,
  phi::FV1, obj::FV2, PI::ProblemInstance,
  Cx::Union{Nothing,Vf},
  Ttilde::Mf, Ttilde_rowmax::Base.RefValue{Ti},
  c::Vector{Vr}, a::Vector{Vf}, b::Vector{Vf}, active::Vector{Bool}, diagjX::Vf,
  debug::Bool=false, BXratio::Tf=1.0
) where {
  Tf<:AbstractFloat,Vf<:AbstractVector{Tf},Mf<:AbstractMatrix{Tf},
  Ti<:Integer,Vi<:AbstractVector{Ti},
  Tr<:Rational{Ti},Vr<:AbstractVector{Tr},
  FV1<:FunctionValue,FV2<:FunctionValue
}
  # initialize
  linsys = 0
  
  # H
  if sum(active) == 0 # case 1
    ctime_H = @elapsed begin
    considx_box!(diagjX, utilde, PI.xlo, PI.xhi)
    t_s = tau/sigma
    if PI.objtype == 2
      @inbounds @simd for i in 1:PI.n
        phi.H[i,i] = obj.H[i,i] + t_s + (sigma * BXratio) * diagjX[i]
      end
    else
      @inbounds @simd for i in 1:PI.n
        phi.H[i,i] = t_s + (sigma * BXratio) * diagjX[i]
      end
    end
    end # ctime_H
    ctime_nd = @elapsed d .= -phi.g ./ diag(phi.H)
  else # case 2
    ctime_H = @elapsed begin
    # prepare linear system
    idx_lo, idx_hi = prepareHphi!(
      xssn, xalm, d,
      ytilde, ybar, ytmp, sig, k0bar, k1bar,
      utilde, ubar,
      sigma, tau,
      phi, obj, PI,
      Ttilde, Ttilde_rowmax, c, a, b, active, diagjX,
      debug, BXratio
    )
    
    # solve linear system
    if PI.objtype == 2
      pobj!(obj, xssn, PI.C, PI.c, Cx)
      diagDinv .= 1.0 ./ (diag(obj.H) .+ tau/sigma .+ (sigma * BXratio) .* diagjX)
    else
      diagDinv .= 1.0 ./ (tau/sigma .+ (sigma * BXratio) .* diagjX)
    end
    diagDinvsqrt .= sqrt.(diagDinv)
    diagDsqrt .= 1.0 ./ diagDinvsqrt
    Dinv_g .= diagDinv .* -phi.g
    @views TT = Ttilde[1:idx_hi,:] #! remember B = A[1:2,:] creates a copy!
    j = zeros(Tf, idx_hi, idx_hi)
    # compute faster (via BLAS.scal!): TT .= TT * Diagonal(diagDinvsqrt)
    @simd for i in 1:PI.n
      @inbounds BLAS.scal!(diagDinvsqrt[i], view(TT, :, i))
    end
    BLAS.gemm!('N', 'T', true, TT, TT, false, j)
    # add σ⁻¹I
    sinv = 1.0 / sigma
    @simd for i in 1:idx_hi
      @inbounds j[i,i] += sinv
    end
    linsys = size(j, 1)
    if debug
      display(sum(j))
    end
    end # ctime_H
    
    ctime_nd = @elapsed begin
    @views LAPACK.potrf!('U', j) #! time it
    diagDsqrt .*= Dinv_g
    sol_M_swm[idx_hi+1:end] .= 0.0
    @views BLAS.gemv!('N', true, TT, diagDsqrt, false, sol_M_swm[1:idx_hi])
    @views LAPACK.potrs!('U', j, sol_M_swm[1:idx_hi]) #! time it
    @views BLAS.gemv!('T', true, TT, sol_M_swm[1:idx_hi], false, sol_n_swm)
    d .= Dinv_g .- (diagDinvsqrt .* sol_n_swm)
    end # ctime_nd
  end # case 2
  return linsys, ctime_H, ctime_nd
end

function ssn!(subproblem_info::Dict,
  xalm::Vf, xssn::Vf, xls::Vf, xhat::Vf, d::Vf,
  Ax::Vector{Vf}, Ad::Vector{Vf},
  ytilde::Vector{Vf}, ybar::Vector{Vf}, sig::Vector{Vi},
  utilde::Vf, ubar::Vf,
  yhat::Vector{Vf}, ytmp::Vector{Vf}, k0bar::Vi, k1bar::Vi,
  uhat::Vf,
  diagDinv::Vf, diagDinvsqrt::Vf, diagDsqrt::Vf, Dinv_g::Vf, sol_M_swm::Vf, sol_n_swm::Vf,
  lambda::Vector{Vf}, mux::Vf, negmux_neg::Vf, negmux_pos::Vf,
  sinvlambda::Vector{Vf}, sinvmux::Vf,
  Atlampmu::Vf, Btlam::Vector{Vf},
  sigma0::Tf, tau0::Tf,
  phi_full::FV1, phi_smw::FV2, phi::Base.RefValue{FunctionValue}, obj::FV3,
  diags::Tuple,
  PI::ProblemInstance,
  Cx::Union{Nothing,Vf}, Ry::Vector{Vf}, Ru::Vf, Rd::Vf,
  Ttilde::Mf, Ttilde_rowmax::Base.RefValue{Ti}, c::Vector{Vr}, a::Vector{Vf}, b::Vector{Vf}, active::Vector{Bool}, diagjX::Vf, mks::Vf,
  maxiter::Ti=100, maxlit::Ti=100, tol::Float64=1e-7,
  ls_method::Symbol=:sw, ls_tau::Tf=0.7, ls_c1::Tf=1e-4, ls_c2::Tf=0.9, stagit::Ti=15, sort_buffer::Tf=0.1,
  lin_solver::Tuple=(:M,1),
  outertol::Tf=1e-8, maxn::Ti=50_000,
  debug::Bool=false, verb::Int64=1, BXratio::Tf=1.0
) where {
  Tf<:AbstractFloat,Vf<:AbstractVector{Tf},Mf<:AbstractMatrix{Tf},
  Ti<:Integer,Vi<:AbstractVector{Ti},
  Tr<:Rational{Ti},Vr<:AbstractVector{Tr},
  FV1<:FunctionValue,FV2<:FunctionValue,FV3<:FunctionValue,
}
  #
  # initialization
  #

  ctime_assign = @elapsed begin
  ssnit::Ti = 0; # newton iterations
  lit::Ti = 0
  lsit::Ti = 0; # linesearch iterations
  sp::Ti = 0; # (intermediate) strict projections value
  sps::Ti = 0; # strict projection count
  psc::Ti = 0; # partial sort count
  linsys::Ti = 0; # lin system size
  absbeta::Ti = 0; # sum(|βl| for l in 1:L)
  sigma = sigma0;
  tau = tau0;
  n::Ti = length(xalm);
  xssn .= xalm;
  stepmin::Tf = ls_tau^(maxlit-3);
  stepmincount::Ti = 0; # stepsize stagnate: small, small, small
  steplast::Tf = 0.0;
  stepcyclecount::Ti = 0; # stepsize cycle: 1.0, small, 1.0, small
  guessquad::Bool = false; # if improve by order of magnitude from previous step
  alpha::Tf = 1.0;
  dCd::Tf = 0.0;
  retcode::Ti = 0;
  phix::Tf = phi[].v;
  gphi_d::Tf = +Inf;
  gobj_d::Tf = +Inf;
  pinfeas, y_Axpb, u_x, dinfeas, relgap, pobj, dobj = diags; # diagnostics
  normgphix = norm(phi[].g); # stopping
  normgphixlast::Tf = Inf; # stopping
  err = +Inf;
  ssnstop::Bool = false;
  re_sort::Bool = true;
  sort_freq::Ti = 1;

  # tracking
  ctime_proj = 0.0; time_proj = 0.0
  ctime_sort = 0.0; time_sort = 0.0
  ctime_vphi = 0.0; time_vphi = 0.0
  time_assign = 0.0
  ctime_gphi = 0.0
  ctime_Hphi = 0.0; time_Hphi = 0.0
  ctime_Ax = 0.0
  ctime_nd = 0.0; time_nd = 0.0

  # stopping
  if tau > 0
    # n small vs m
    resnorm_buffer_length = 5
    max_cycle_length = 5

    # n large vs m
    resnorm_buffer_length = 10
    max_cycle_length = 8
  else
    resnorm_buffer_length = 10
    max_cycle_length = 8
  end
  resnorm_buffer = zeros(resnorm_buffer_length)
  resnorm_ratio = zeros(resnorm_buffer_length)
  stagit_check = SP.stagit
  cycleit_check = SP.stagit
  cycletol = 1e-4
  end # ctime_assign

  #
  # printing
  #

  absbeta = sum(k1bar .- k0bar)
  subproblem_info[:ctime_diagnostics] += @elapsed begin
  if verb >= 2
    print(SUB_line(0, pinfeas, y_Axpb, u_x, dinfeas, pobj, relgap, phix, normgphix, linsys, absbeta, 0.0, 0, active, 0.0, true, false))
  end; end # ctime_diagnostics
  dinfeaslast = +Inf;
  # for l in 1:PI.L
  #   println("k0    k1    ($l) = ", get_k0k1!(ytilde[l][sig[l]], PI.k[l]))
  #   println("k0bar k1bar ($l) = ", get_k0k1!(ybar[l][sig[l]], PI.k[l]))
  # end
  
  #
  # perform ssn iters
  #
  
  ssn_start = time();
  ssnit_start = Inf;
  ssnit_stop = Inf;
  
  while true
    ssnit_start = time();
    ssnit += 1;
    steplast = alpha;
    normgphixlast = normgphix;

    if debug
      println("========== ssn ==========")
      println("    dobj = $(phi[].d)")
      println("    pobj = $(phi[].v)")
      println("    lam  = $lambda")
      println("    b  = $(PI.b)")
      println("    b*lam = $(lambda' * PI.b)")
    end

    #
    # get newton direction
    #

    # if lin_solver[1] == :M && lin_solver[2] == 1
    if n <= maxn && PI.n <= 10sum(PI.m)*sum(PI.k)/sum(PI.n) # either smw or full
      if sum(active) == 0 && n <= 2sum(PI.m)
      # if true
        # println("using full")
        if debug; println("using full"); end
        time_Hphi, time_nd = nd!(
          xssn, xalm, d,
          ytilde, ybar, ytmp, yhat, sig, k0bar, k1bar,
          utilde, ubar, uhat,
          sigma, tau,
          phi_full, obj, PI,
          Cx,
          Ttilde, Ttilde_rowmax, c, a, b, active, diagjX,
          debug, BXratio
        )
        ctime_assign += @elapsed begin
        linsys = PI.n
        phi[] = phi_full
        phi_smw.v = phi_full.v
        phi_smw.g .= phi_full.g
        phi_smw.d = phi_full.d
        end # ctime_assign
        subproblem_info[:avg_nd_size] += linsys
        subproblem_info[:absbeta] += -sum(k0bar[active].-k1bar[active])/sum(active)
        subproblem_info[:absalpha] += sum(k0bar[active])/sum(active)
      elseif sum(absbeta) <= n
        # println("using smw")
        if debug; println("using smw"); end
        linsys, time_Hphi, time_nd = nd_swm!(
          xssn, xalm, d,
          ytilde, ybar, ytmp, yhat, sig, k0bar, k1bar,
          utilde, ubar, uhat,
          diagDinv, diagDinvsqrt, diagDsqrt, Dinv_g, sol_M_swm, sol_n_swm,
          sigma, tau,
          phi_smw, obj, PI,
          Cx,
          Ttilde, Ttilde_rowmax, c, a, b, active, diagjX,
          debug, BXratio
        )
        ctime_assign += @elapsed begin
        phi[] = phi_smw
        phi_full.v = phi_smw.v
        phi_full.g .= phi_smw.g
        phi_full.d = phi_smw.d
        end # ctime_assign
        subproblem_info[:avg_nd_size] += linsys
        subproblem_info[:absbeta] += -sum(k0bar[active].-k1bar[active])/sum(active)
        subproblem_info[:absalpha] += sum(k0bar[active])/sum(active)
      # elseif lin_solver[1] == :n && lin_solver[2] == 1
      elseif sum(absbeta) > n
        # println("using full")
        if debug; println("using full"); end
        time_Hphi, time_nd = nd!(
          xssn, xalm, d,
          ytilde, ybar, ytmp, yhat, sig, k0bar, k1bar,
          utilde, ubar, uhat,
          sigma, tau,
          phi_full, obj, PI,
          Cx,
          Ttilde, Ttilde_rowmax, c, a, b, active, diagjX,
          debug, BXratio
        )
        ctime_assign += @elapsed begin
        linsys = PI.n
        phi[] = phi_full
        phi_smw.v = phi_full.v
        phi_smw.g .= phi_full.g
        phi_smw.d = phi_full.d
        end # ctime_assign
        subproblem_info[:avg_nd_size] += linsys
        subproblem_info[:absbeta] += -sum(k0bar[active].-k1bar[active])/sum(active)
        subproblem_info[:absalpha] += sum(k0bar[active])/sum(active)
      else
        throw()
      end
    else # too big for full
      # println("too big for full: using swm")
      linsys, time_Hphi, time_nd = nd_swm!(
        xssn, xalm, d,
        ytilde, ybar, ytmp, yhat, sig, k0bar, k1bar,
        utilde, ubar, uhat,
        diagDinv, diagDinvsqrt, diagDsqrt, Dinv_g, sol_M_swm, sol_n_swm,
        sigma, tau,
        phi_smw, obj, PI,
        Cx,
        Ttilde, Ttilde_rowmax, c, a, b, active, diagjX,
        debug, BXratio
      )
      ctime_assign += @elapsed begin
      phi[] = phi_smw
      phi_full.v = phi_smw.v
      phi_full.g .= phi_smw.g
      phi_full.d = phi_smw.d
      end # ctime_assign
      subproblem_info[:avg_nd_size] += linsys
      subproblem_info[:absbeta] += -sum(k0bar[active].-k1bar[active])/sum(active)
      subproblem_info[:absalpha] += sum(k0bar[active])/sum(active)
    end
    ctime_Hphi += time_Hphi
    ctime_nd += time_nd

    #
    # update Ad, ⟨∇φ,d⟩, ⟨∇f,d⟩, dCd
    #
    
    ctime_Ax += @elapsed begin
    for l in eachindex(PI.m) # Threads.@threads
      # Ad[l] .= PI.A[l] * d
      # gemv!(tA, alpha, A, x, beta, y) => y = alpha*A*x + beta*y
      @views BLAS.gemv!('N', true, PI.A[l], d, false, Ad[l])
    end; end # ctime_Ax
    gphi_d = dot(phi[].g, d)
    gobj_d = dot(obj.g, d)
    if PI.objtype == 2
      dCd = d' * PI.C * d
    end

    #
    # linesearch
    #

    fudge = 1e-8
    alpha = 1.0
    if mod(ssnit, sort_freq) == 0
      re_sort = true
    else
      re_sort = false
    end
    if ls_method == :sw
      alpha, lit, phix, gphix, sp, time_proj, time_sort, time_vphi, time_assign, psc, lsretcode = strongwolfe!(
        xls, xssn, xalm, xhat, d,
        Ax, Ad,
        ytilde, ybar, yhat, ytmp, sig, k0bar, k1bar, active, mks,
        sinvlambda, sinvmux,
        utilde, ubar, uhat,
        phix, gphi_d,
        dCd, gobj_d,
        phi[], obj,
        PI,
        Cx,
        sigma, tau,
        alpha, ls_tau, ls_c1, ls_c2, maxlit, re_sort, sort_buffer,
        fudge,
        debug, BXratio
      )
    elseif ls_method == :bt
      alpha, lit, phix, sp, time_proj, time_sort, time_vphi, time_assign, psc = backtracking!(
        xls, xssn, xalm, d,
        Ax, Ad,
        ytilde, ybar, yhat, ytmp, sig, k0bar, k1bar, active, mks,
        sinvlambda, sinvmux,
        utilde, ubar, uhat,
        phix, gphi_d,
        phi[], obj,
        PI,
        Cx,
        sigma, tau,
        alpha, ls_tau, ls_c1, maxlit, re_sort, sort_buffer,
        fudge,
        debug, BXratio
      )
    elseif ls_method == :ww
      alpha, lit, phix, gphix, sp, time_proj, time_sort, time_vphi, time_assign, psc = weakwolfe!(
        xls, xssn, xalm, xhat, d,
        Ax, Ad,
        ytilde, ybar, yhat, ytmp, sig, k0bar, k1bar, active, mks,
        sinvlambda, sinvmux,
        utilde, ubar, uhat,
        phix, gphi_d,
        dCd, gobj_d,
        phi[], obj,
        PI,
        Cx,
        sigma, tau,
        alpha, ls_tau, ls_c1, maxlit, re_sort, sort_buffer,
        fudge,
        debug, BXratio
      )
    end
    lsit += lit
    sps += sp
    ctime_proj += time_proj
    ctime_sort += time_sort
    ctime_vphi += time_vphi
    ctime_assign += time_assign
    subproblem_info[:avg_step_size] += alpha
    subproblem_info[:partial_sort_count] += psc
    
    #
    # updates
    #

    # update Ax
    # note: y(x) and u(x) are updated in ls
    ctime_Ax += @elapsed begin
    for l in 1:PI.L # Threads.@threads 
      # Ax[l] .= PI.A[l] * x
      # gemv!(tA, alpha, A, x, beta, y) => y = alpha*A*x + beta*y
      # @views BLAS.gemv!('N', true, PI.A[l], d, false, Ad[l])
      Ax[l] .+= alpha .* Ad[l] #! riskier? but seems to work
    end; end # time_Ax
    
    # update gradient and stopping value #! computes A' * yhat
    ctime_gphi += @elapsed gphi!(phi[], obj, xssn, xalm, yhat, uhat, sigma, tau, PI, Cx, BXratio)
    normgphix = norm(phi[].g)
    resnorm_buffer[mod1(ssnit,resnorm_buffer_length)] = normgphix
    
    # diagnostics
    dinfeaslast = dinfeas
    subproblem_info[:ctime_diagnostics] += @elapsed begin
    (pinfeas, y_Axpb, u_x, dinfeas, pobj, relgap, phix, normgphix) = calc_SUB_line(
      xssn, xalm, ybar, ubar, lambda, mux, negmux_neg, negmux_pos, Ax, Btlam, Atlampmu, phi[], obj, sigma, tau, PI,
      Cx, Ry, Ru, Rd, sig,
    )
    end # ctime_diagnostics
    err = max(abs(pinfeas), abs(dinfeas), abs(relgap))
    
    #
    # end iteration
    #
    
    ssnit_stop = time()
    

    #
    # stopping
    #

    # solve alm problem
    if err <= outertol
      retcode = 1
      ssnstop = true
    end

    # solve subproblem
    if normgphix <= tol
      retcode = 2
      ssnstop = true
    end

    # iteration limit
    if ssnit >= maxiter
      retcode = 3
      ssnstop = true
    end

    # stagnation
    if (ssnit >= resnorm_buffer_length) && (ssnit >= stagit_check)
      resnorm_ratio .= resnorm_buffer[mod1(ssnit,resnorm_buffer_length)]./resnorm_buffer
      # println(resnorm_ratio)
      # resnorm
      if tau > 0
        if PI.n >= 100_000
          if (
            (minimum(resnorm_ratio) > 1/2) &&
            (maximum(resnorm_ratio) < 2)
          )
            println("stopped phi stagnation!")
            if debug; println("stopped stagnation!"); end
            ssnstop = true
            retcode = 4
          end
        end
      else
        if (
          (minimum(resnorm_ratio) > 1/2) &&
          (maximum(resnorm_ratio) < 2)
        )
          println("stopped phi stagnation!")
          if debug; println("stopped stagnation!"); end
          ssnstop = true
          retcode = 4
        end
      end
    end

    # cycling
    if (ssnit > resnorm_buffer_length) && (ssnit > cycleit_check)
      # num = resnorm_buffer[mod1(ssnit,resnorm_buffer_length)]
      # den = resnorm_buffer[mod1.(ssnit-max_cycle_length-1:ssnit-2,resnorm_buffer_length)]
      # v = abs.(1 .- (num ./ den))
      if debug
        println("checking buffer")
        println("num = $num")
        println("den = $den")
        println("v   = $v")
      end
      # if minimum(v) <= cycletol
      if minimum(
        abs.( 1.0 .- 
          resnorm_buffer[mod1(ssnit,resnorm_buffer_length)]
          ./
          resnorm_buffer[mod1.(ssnit-max_cycle_length-1:ssnit-2,resnorm_buffer_length)]
        )
      ) <= cycletol
        println("stopped cycling!")
        if debug
          println("stop cycling!")
        end
        ssnstop = true
        retcode = 6
      end
    end

    if time()-subproblem_info[:start_time] > AP.maxtime
      ssnstop = true
      retcode = 7
    end

    # stop
    if ssnstop
      absbeta = sum(k1bar .- k0bar)
      subproblem_info[:ctime_diagnostics] += @elapsed begin
      if verb >= 2
        print(SUB_line(ssnit, pinfeas, y_Axpb, u_x, dinfeas, pobj, relgap, phix, normgphix, linsys, absbeta, alpha, lit, active, ssnit_stop-ssnit_start, false, true))
      end; end # ctime_diagnostics
      if retcode == 1 || retcode == 2
        subproblem_info[:ctime_diagnostics] += @elapsed begin
        if verb >= 2; print(SUB_line_term(dinfeas, normgphix, tol, true)); end
        end # ctime_diagnostics
      else
        subproblem_info[:ctime_diagnostics] += @elapsed begin
        if verb >= 2; print(SUB_line_term(dinfeas, normgphix, tol, false)); end
        end # ctime_diagnostics
      end
      break
    end
    
    #
    # printing
    #

    absbeta = sum(k1bar .- k0bar)
    subproblem_info[:ctime_diagnostics] += @elapsed begin
    if verb >= 2
      print(SUB_line(ssnit, pinfeas, y_Axpb, u_x, dinfeas, pobj, relgap, phix, normgphix, linsys, absbeta, alpha, lit, active, ssnit_stop-ssnit_start, false, false))
    end; end # ctime_diagnostics
  end # while true

  # stop timer
  ssn_stop = time()
  
  # subproblem tracking
  ctime_assign += @elapsed begin
  subproblem_info[:iter] = ssnit
  subproblem_info[:lsiter] = lsit
  subproblem_info[:avg_step_size] /= ssnit
  subproblem_info[:sp] = sps
  subproblem_info[:err] = err
  subproblem_info[:status] = retcode
  subproblem_info[:time] = ssn_stop-ssn_start
  subproblem_info[:ctime_proj] = ctime_proj
  subproblem_info[:ctime_sort] = ctime_sort
  subproblem_info[:ctime_nd] = ctime_nd
  subproblem_info[:ctime_Hphi] = ctime_Hphi
  subproblem_info[:ctime_gphi] = ctime_gphi
  subproblem_info[:ctime_vphi] = ctime_vphi
  subproblem_info[:ctime_Ax] = ctime_Ax
  subproblem_info[:avg_nd_size] = subproblem_info[:avg_nd_size] / ssnit
  subproblem_info[:absbeta] = subproblem_info[:absbeta] / ssnit
  subproblem_info[:absalpha] = subproblem_info[:absalpha] / ssnit
  subproblem_info[:normg] = normgphix
  subproblem_info[:pinfeas] = pinfeas
  subproblem_info[:dinfeas] = dinfeas
  subproblem_info[:relgap] = relgap
  subproblem_info[:pobj] = pobj
  subproblem_info[:dobj] = dobj
  end # ctime_assign
  subproblem_info[:ctime_assign] += ctime_assign
end
