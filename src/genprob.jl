function ObjectiveInstance(seed::Ti, objtype::Ti, n::Ti) where Ti<:Integer
  rng = Xoshiro(seed)
  c = 1randn(rng, n)
  if objtype == 2
    Cdiag = abs.(1randn(rng, n))
    # Cdiag = ones(n)
    C = Diagonal(Cdiag)
    Cinv = Diagonal(1.0 ./ Cdiag)
  else
    C = nothing
    Cinv = nothing
  end
  return ObjectiveInstance(seed, objtype, C, Cinv, c)
end
function ConstraintInstance(seed::Ti, n::Ti, m::Vi, k::Vi, scale_row::Bool, scale_col::Bool, t::Tf) where {Tf<:AbstractFloat,Ti<:Integer,Vi<:AbstractVector{Ti}}
  rng = Xoshiro(seed)
  @assert(scale_row + scale_col <= 1)
  # @assert(t <= 1.0 && t >= 0.0)
  A = Matrix{Tf}[]
  b = Vector{Tf}[]
  L = length(m)
  for l in eachindex(m)
    # A matrix
    Al = 10randn(rng, m[l], n)
    # Al = rand(rng, -5:5, m[l], n)./1.0
    if scale_row
      for i in 1:m[l]
        # Al[i,:] ./= norm(Al[i,:], 2)
        Al[i,:] ./= norm(Al[i,:], Inf)
      end
    end
    if scale_col
      for j in 1:n
        # Al[:,j] ./= norm(Al[:,j], 2)
        Al[:,j] ./= norm(Al[:,j], Inf)
      end
    end

    # b vector
    if scale_col
      # bl = sqrt(n) .* randn(m[l])
      bl = randn(rng, m[l])
    end
    if scale_row
      # bl = sqrt(m[l]) .* randn(m[l])
      bl = randn(rng, m[l])
    end
    # bl = rand(-5:5,m[l])./1.0
    push!(A, Al)
    push!(b, bl)
  end
  
  # box constraint
  xlo = fill(-1.0, n);
  xhi = fill(+1.0, n);
  # xlo = fill(-1.0, n)./10.0;
  # xhi = fill(+1.0, n)./10.0;
  # xlo = -rand(n);
  # xhi = +rand(n);
  xloidx = fill(true, n);
  xhiidx = fill(true, n);
  Xunconstr = false;
  
  # max-k-sum rhs
  # nfeas = Int(ceil(sqrt(n)))
  nfeas = Int(ceil(log(n+sum(m))))
  r_hard = [zeros(m[l], nfeas) for l in eachindex(m)]
  r_wtd = [zeros(m[l], nfeas) for l in eachindex(m)]
  for j in 1:nfeas
    # xfeas = rand(rng, hcat(xlo, xhi), n) # on boundary
    xfeas = xlo .+ (xhi .- xlo) .* rand(rng, n)
    for l in eachindex(m)
      Al = A[l]
      bl = b[l]
      r_hard[l][:,j] .= Al * xfeas .+ bl
    end
  end
  sort!.(r_hard, rev=true, dims=1)
  mks_hard = [sum(r_hard[l][1:k[l],:], dims=1) for l in eachindex(m)]
  mks_hard = vcat(mks_hard...)
  mks_hard_min = mks_hard[:,argmin(vec(sum(mks_hard,dims=1)))]
  mks_hard_max = mks_hard[:,argmax(vec(sum(mks_hard,dims=1)))]
  mks_wtd = t .* mks_hard_max .+ (1-t) .* mks_hard_min
  println("                      t: $t")
  println("sum_{i=1:k}(rl[i]) easy: $(mks_hard_max)")
  println("sum_{i=1:k}(rl[i]) hard: $(mks_hard_min)")
  println("sum_{i=1:k}(rl[i]) used: $(mks_wtd)")
  
  # B
  Bt = Function[]
  for l in eachindex(m)
    Btl = Atmatinv_action(m[l], k[l])
    push!(Bt, Btl)
  end

  # CI
  return ConstraintInstance(
    seed, n, m, L, A, b, k, mks_wtd, Bt, xlo, xhi, xloidx, xhiidx, Xunconstr,
  )
end

function make_quantile_regression(PI::ProblemInstance, epi::Vector{Bool})
  ## setup
  Tf = eltype(PI.c)
  L = PI.L
  m = PI.m
  nx = PI.n
  nt = sum(epi)

  ## epigraph reformulation
  # constraint: T_k(y) <= t <==> Σ₁ᵏ y_(i) - t/k <= 0
  DT = eltype(PI.A)
  A = DT[]
  DT = eltype(PI.r)
  r = DT[]
  for l in 1:L
    if epi[l]
      push!(A, hcat(PI.A[l], -ones(Tf, m[l]))) # ok
      push!(r, 0.0)
    else
      push!(A, PI.A[l])
      push!(r, PI.r[l])
    end
  end
  xlo = [PI.xlo; fill(-Inf, nt)]
  xhi = [PI.xhi; fill(+Inf, nt)]
  xloidx = [PI.xloidx; fill(false, nt)]
  xhiidx = [PI.xhiidx; fill(false, nt)]
  # variables
  n = nx + nt
  # objective
  c = [PI.c; ones(Tf, nt)]
  if !isnothing(PI.C)
    if isa(PI.C, Diagonal)
      C = Diagonal([diag(C); zeros(Tf, nt)])
      Cinv = Diagonal([1.0 ./ diag(C); zeros(Tf, nt)])
    end
  else
    C = nothing
    Cinv = nothing
  end

  return ProblemInstance(
    PI.conseed, PI.objseed, n, PI.m, PI.L,
    A, PI.b, PI.k, r, PI.Bt,
    xlo, xhi, xloidx, xhiidx, PI.Xunconstr,
    PI.objtype, C, Cinv, c, PI.scale
  )
end
