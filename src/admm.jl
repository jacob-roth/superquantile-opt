"""
3step admm
"""
function admm3(
  x0::Vf, lambda0::Vector{Vf}, sigma0::Tf, tau0::Tf, PI::ProblemInstance,
  obj::FunctionValue,
  ditmax::Integer=10000, tol::Tf=5e-3, largedualstep::Bool=false, debug::Bool=false,
) where {Tf<:AbstractFloat, Vf<:AbstractVector{Tf}}
  # initialize
  x = zeros(Tf, PI.n)
  xs = zeros(Tf, PI.n, ditmax+1)
  xavg = zeros(Tf, PI.n)
  xprev = zeros(Tf, PI.n)
  y = [zeros(Tf, PI.m[l]) for l in 1:PI.L]
  yavg = [zeros(Tf, PI.m[l]) for l in 1:PI.L]
  yy = [zeros(Tf, PI.m[l]) for l in 1:PI.L]
  yhat = [zeros(Tf, PI.m[l]) for l in 1:PI.L]
  ytmp = [zeros(Tf, PI.m[l]) for l in 1:PI.L]
  z = zeros(Tf, PI.n)
  zavg = zeros(Tf, PI.n)
  Ax = [zeros(Tf, PI.m[l]) for l in 1:PI.L]
  lambda = [zeros(Tf, PI.m[l]) for l in 1:PI.L]
  mu = zeros(Tf, PI.n)
  negmux_neg = zeros(Tf, PI.n); # [-mux]^-
  negmux_pos = zeros(Tf, PI.n); # [-mux]^+
  xlo = zeros(Tf, PI.n)
  xhi = zeros(Tf, PI.n)
  xlo .= PI.xlo
  xhi .= PI.xhi
  x .= x0
  xprev .= x0
  z .= x0
  rhs = zeros(Tf, PI.n)
  b_M = zeros(Tf, sum(PI.m))

  # helper: matvec
  Ax = [zeros(Tf, m[l]) for l in eachindex(m)];
  Ad = [zeros(Tf, m[l]) for l in eachindex(m)];
  Atlampmu = zeros(Tf, n); # A' * lambda + mu
  Btlam = [zeros(Tf, m[l]) for l in eachindex(m)]; # Atildeinv' * lambda
  Cx = zeros(Tf, n);
  
  
  dit = 1
  pwin = 0
  dwin = 0
  id = 1
  pinfeas, dinfeas, relgap, pobj, dobj_, = Inf, Inf, Inf, Inf, Inf
  sigma = sigma0
  for l in 1:PI.L
    lambda[l] .= lambda0[l]
    Ax[l] .= PI.A[l] * x
    yy[l] = Ax[l] .+ PI.b[l] .- lambda[l] ./ sigma
    sig_yy = sortperm(yy[l], rev=true)
    # @views project_maxksum!(y[l][sig_yy], yy[l][sig_yy], PI.r[l], PI.k[l])
    @views _, (k0l, k1l), _ = @views project_maxksum_snake!(y[l][sig_yy], yy[l][sig_yy], PI.r[l], PI.k[l], sum(yy[l][sig_yy[1:PI.k[l]]]) > PI.r[l], false) #! here ybar[l] is assumed ("incorrectly") to be sorted
  end
  A = vcat(PI.A...)
  b = vcat(PI.b...)
  C = PI.C
  c = PI.c
  
  #! only when τ = 0
  if sum(PI.m) < PI.n
    if PI.objtype == 1
      println("M < n & linear obj ==> factor (I + AAᵀ)")
      AAt = A * A'
      I_AAt = AAt + Diagonal(ones(sum(PI.m)))
      f_I_AAt = LAPACK.potrf!('U', I_AAt)
    elseif PI.objtype == 2
      println("M < n & quadratic (diagonal) obj ==> factor (σ⁻¹I + A D⁻¹ Aᵀ)")
      D = C + Diagonal(sigma .* ones(PI.n))
      Dinv = Diagonal(1.0 ./ diag(D))
      sinvI_ADAt = Diagonal((1.0/sigma) .* ones(sum(PI.m))) + A * Dinv * A'
      f_sinvI_ADAt = LAPACK.potrf!('U', sinvI_ADAt)
    end
  else
    AtA = A'*A
    if PI.objtype == 1
      println("M > n & linear obj ==> factor (I + AᵀA)")
      I_AtA = AtA .+ Diagonal(sigma .* ones(PI.n))
      f_I_AtA = LAPACK.potrf!('U', I_AtA)
    elseif PI.objtype == 2
      println("M > n & quadratic (diagonal) obj ==> factor (σ⁻¹C + I + AᵀA)")
      sC_I_AtA = C ./ sigma .+ (Diagonal(ones(PI.n)) .+ AtA)
      f_sC_I_AtA = LAPACK.potrf!('U', sC_I_AtA)
    end
  end
  
  while dit <= ditmax
    t_start = time()

    # x update: O(n^2 + nm) per step
    rhs .= -(c .- A' * vcat(lambda...) .- sigma .* A' * (-b + vcat(y...)) .- sigma .* z .- mu)
    # rhs .= -(c .- A' * lambda[1] .- sigma .* A' * (-b + y[1]) .- sigma .* z .- mu)
    if PI.objtype == 1
      if sum(PI.m) < PI.n
        # "x = rhs - A' * ((I + AAt) \ (A * rhs))"
        b_M .= A * rhs
        LAPACK.potrs!('U', f_I_AAt[1], b_M)
        x .= rhs .- A' * b_M
      else
        # "x = sigma .* ((I + A'*A) \ (sigma .* rhs))"
        x .= sigma .* rhs
        x .= LAPACK.potrs!('U', f_I_AtA[1], x)
        x .*= sigma
      end
    else
      if sum(PI.m) < PI.n
        # "x = Dinv * rhs - Dinv A' * ((sinvI + A*Dinv*A') \ (A * (Dinv * rhs)))"
        x .= Dinv * rhs
        b_M .= A * x
        LAPACK.potrs!('U', f_sinvI_ADAt[1], b_M)
        x .= Dinv * rhs .- Dinv * (A' * b_M)
      else
        # "x = ((sigma⁻¹C + I + A'*A) \ (rhs ./ sigma))"
        x .= rhs ./ sigma
        x .= LAPACK.potrs!('U', f_sC_I_AtA[1], x)
      end
    end

    # y and lambda update: O(mn)
    for l in 1:PI.L
      Ax[l] .= PI.A[l] * x # O(mn)
      yy[l] = Ax[l] .+ PI.b[l] .- lambda[l] ./ sigma
      sig_yy = sortperm(yy[l], rev=true)
      @views _, (k0l, k1l), _ = @views project_maxksum_snake!(y[l][sig_yy], yy[l][sig_yy], PI.r[l], PI.k[l], sum(yy[l][sig_yy[1:PI.k[l]]]) > PI.r[l], false) #! here ybar[l] is assumed ("incorrectly") to be sorted
      if largedualstep
        lambda[l] .= min.(lambda[l] .+ sqrt(ℯ)sigma .* (y[l] .- Ax[l] .- PI.b[l]), 0.0)
      else
        lambda[l] .= min.(lambda[l] .+ sigma .* (y[l] .- Ax[l] .- PI.b[l]), 0.0)
      end
    end

    # z and mu update # O(n)
    @views project_box!(z, x .- mu ./ sigma, xlo, xhi)
    if largedualstep
      mu .+= sqrt(ℯ)sigma .* (z .- x)
    else
      mu .+= sigma .* (z .- x)
    end

    # avg
    xavg .= ((dit-1) .* xavg .+ x) ./ dit
    for l in 1:PI.L
      yavg[l] .= ((dit-1) .* yavg[l] .+ y[l]) ./ dit
    end
    zavg .= ((dit-1) .* zavg .+ x) ./ dit

    # diagnostics and stopping
    if mod(dit, 25) == 0
      (pinfeas_avg, dinfeas_avg, relgap_avg, pobj_avg, dobj_avg, _) = calc_ADMM_line(xavg, xavg, yavg, zavg, lambda, mu, negmux_neg, negmux_pos, Ax, Btlam, Atlampmu, obj, sigma, 0.0, PI)
      (pinfeas_, dinfeas_, relgap_, pobj_, dobj_, _) = calc_ADMM_line(x, x, y, z, lambda, mu, negmux_neg, negmux_pos, Ax, Btlam, Atlampmu, obj, sigma, 0.0, PI)
      if max(pinfeas_avg, dinfeas_avg, relgap_avg) <= max(pinfeas_, dinfeas_, relgap_)
        pinfeas, dinfeas, relgap, objp, objd = pinfeas_avg, dinfeas_avg, relgap_avg, pobj_avg, dobj_avg
        id = 1
      else
        pinfeas, dinfeas, relgap, objp, objd = pinfeas_, dinfeas_, relgap_, pobj_, dobj_
        id = 2
      end
      t_stop = time()
      if dit == 1 || mod(dit, 25) == 0
        print(ADMM_line(dit, pinfeas, dinfeas, relgap, objp, objd, sigma, 0.0, id, t_stop-t_start, true))
      else
        print(ADMM_line(dit, pinfeas, dinfeas, relgap, objp, objd, sigma, 0.0, id, t_stop-t_start, false))
      end
    end
    xs[:,dit] .= x
    dit += 1

    if max(abs(pinfeas), abs(dinfeas), abs(relgap)) <= tol
      break
    end
  
  end
  if id == 1
    (xbest, ybest, zbest) = (xavg, yavg, zavg)
  elseif id == 2
    (xbest, ybest, zbest) = (x, y, z)
  end
  res = Dict()
  res[:xbest] = xbest
  res[:ybest] = ybest
  res[:zbest] = zbest
  res[:lambda] = lambda
  res[:mu] = mu
  res[:pinfeas] = pinfeas
  res[:dinfeas] = dinfeas
  res[:relgap] = relgap
  res[:pobj] = pobj
  res[:dobj_] = dobj_
  res[:xs] = xs[:,1:dit]
  return res
end

"""
majorize admm
"""
function admmM(
  x0::Vf, lambda0::Vector{Vf}, sigma0::Tf, tau0::Tf, PI::ProblemInstance,
  obj::FunctionValue,
  ditmax::Integer=10000, tol::Tf=5e-3, largedualstep::Bool=false, debug::Bool=false,
) where {Tf<:AbstractFloat, Vf<:AbstractVector{Tf}}
  # initialize
  x = zeros(Tf, PI.n)
  xs = zeros(Tf, PI.n, ditmax+1)
  xavg = zeros(Tf, PI.n)
  xprev = zeros(Tf, PI.n)
  y = [zeros(Tf, PI.m[l]) for l in 1:PI.L]
  yavg = [zeros(Tf, PI.m[l]) for l in 1:PI.L]
  yy = [zeros(Tf, PI.m[l]) for l in 1:PI.L]
  z = zeros(Tf, PI.n)
  zavg = zeros(Tf, PI.n)
  Ax = [zeros(Tf, PI.m[l]) for l in 1:PI.L]
  lambda = [zeros(Tf, PI.m[l]) for l in 1:PI.L]
  mu = zeros(Tf, PI.n)
  xlo = zeros(Tf, PI.n)
  xhi = zeros(Tf, PI.n)
  xlo .= PI.xlo
  xhi .= PI.xhi
  x .= x0
  xprev .= x0
  z .= x0
  rhs = zeros(Tf, PI.n)
  b_M = zeros(Tf, sum(PI.m))
  dit = 1
  pwin = 0
  dwin = 0
  id = 1
  pinfeas, dinfeas, relgap, pobj, dobj_, = Inf, Inf, Inf, Inf, Inf
  sigma = sigma0
  for l in 1:PI.L
    lambda[l] .= lambda0[l]
    Ax[l] .= PI.A[l] * x
    yy[l] = Ax[l] .+ PI.b[l] .- lambda[l] ./ sigma
    sig_yy = sortperm(yy[l], rev=true)
    @views project_maxksum!(y[l][sig_yy], yy[l][sig_yy], PI.r[l], PI.k[l])
  end
  A = vcat(PI.A...)
  b = vcat(PI.b...)
  
  eigmaxAtA = 0.0
  if sum(PI.m) < PI.n
    _,S,_ = svd(vcat(PI.A...))
    eigmaxAtA = max(eigmaxAtA, maximum(S)^2)
    #! cannot compute singular values of submatrices
  end
  Xnone = all(isinf.(PI.xlo)) && all(isinf.(PI.xhi))
  if Xnone
    Q = Diagonal(sigma*eigmaxAtA .* ones(n))
  else
    Q = Diagonal(sigma .* (ones(n) .+ eigmaxAtA))
  end
  
  while dit <= ditmax
    t_start = time()

    # x update: O(n + nm)
    if PI.objtype == 1
    else
      if Xnone
        rhs .= -(c .- A' * vcat(lambda...) .- sigma .* A' * (-b + vcat(y...)))
      else
        rhs .= -(c .- A' * vcat(lambda...) .- sigma .* A' * (-b + vcat(y...)) .- sigma .* z .- mu)
      end
      x .= (1.0 ./ (diag(Q) .+ diag(C))) .* rhs
    end

    # y and lambda update: O(mn)
    for l in 1:PI.L
      Ax[l] .= PI.A[l] * x # O(mn)
      yy[l] = Ax[l] .+ PI.b[l] .- lambda[l] ./ sigma
      sig_yy = sortperm(yy[l], rev=true)
      @views project_maxksum!(y[l][sig_yy], yy[l][sig_yy], PI.r[l], PI.k[l])
      if largedualstep
        lambda[l] .+= sqrt(ℯ)sigma .* (y[l] .- Ax[l] .- PI.b[l])
      else
        lambda[l] .+= sigma .* (y[l] .- Ax[l] .- PI.b[l])
      end
    end

    # z and mu update: O(n)
    @views project_box!(z, x .- mu ./ sigma, xlo, xhi)
    if largedualstep
      mu .+= sqrt(ℯ)sigma .* (z .- x)
    else
      mu .+= sigma .* (z .- x)
    end

    # avg
    xavg .= ((dit-1) .* xavg .+ x) ./ dit
    for l in 1:PI.L
      yavg[l] .= ((dit-1) .* yavg[l] .+ y[l]) ./ dit
    end
    zavg .= ((dit-1) .* zavg .+ x) ./ dit

    # diagnostics and stopping
    if mod(dit, 25) == 0
      
      (pinfeas_avg, dinfeas_avg, relgap_avg, pobj_avg, dobj_avg, _) = calc_ADMM_line(xavg, xavg, yavg, zavg, lambda, mu, negmux_neg, negmux_pos, Ax, Btlam, Atlampmu, obj, sigma, 0.0, PI)
      (pinfeas_, dinfeas_, relgap_, pobj_, dobj_, _) = calc_ADMM_line(x, x, y, z, lambda, mu, negmux_neg, negmux_pos, Ax, Btlam, Atlampmu, obj, sigma, 0.0, PI)
      if max(pinfeas_avg, dinfeas_avg, relgap_avg) <= max(pinfeas_, dinfeas_, relgap_)
        pinfeas, dinfeas, relgap, objp, objd = pinfeas_avg, dinfeas_avg, relgap_avg, pobj_avg, dobj_avg
        id = 1
      else
        pinfeas, dinfeas, relgap, objp, objd = pinfeas_, dinfeas_, relgap_, pobj_, dobj_
        id = 2
      end
      t_stop = time()
      if dit == 1 || mod(dit, 25) == 0
        print(ADMM_line(dit, pinfeas, dinfeas, relgap, objp, objd, sigma, 0.0, id, t_stop-t_start, true))
      else
        print(ADMM_line(dit, pinfeas, dinfeas, relgap, objp, objd, sigma, 0.0, id, t_stop-t_start, false))
      end
    end
    xs[:,dit] .= x
    dit += 1

    if max(abs(pinfeas), abs(dinfeas), abs(relgap)) <= tol
      break
    end
  
  end
  if id == 1
    (xbest, ybest, zbest) = (xavg, yavg, zavg)
  elseif id == 2
    (xbest, ybest, zbest) = (x, y, z)
  end
  return (xbest, ybest, zbest), (lambda, mu), (pinfeas, dinfeas, relgap, pobj, dobj_), xs[:,1:dit]
end
