#! should rename to "pfp" == primal feasibility pivoting if it turns out this is true
function project_maxksum_snake!(
  xbarsort::AbstractVector{Tfr}, x0sort::AbstractVector{Tfr}, r::Tfr, k::Ti,
  active::Bool, x0prepop::Bool=false, verb::Bool=false, hist::Bool=false,
) where {Tfr<:Union{AbstractFloat,Rational},Ti<:Integer}
time_init_1 = @elapsed begin
  # inactive
  n = length(x0sort)
  if !active
    case = 0
    if !x0prepop
      @simd for i in 1:n
        @inbounds xbarsort[i] = x0sort[i]
      end
    end
    if hist
      return 0, get_k0k1(xbarsort, k), 0, (0.0, 0.0, 0.0), [(0,0)], case, 0.0
    else
      return 0, get_k0k1(xbarsort, k), 0, (0.0, 0.0, 0.0), case, 0.0
    end
  end

  # initialize
  solved::Bool = false
  n = length(x0sort)
  s0::Tfr = sum(view(x0sort, 1:k))
  lam::Tfr = 0
  case::Ti = 0
  # println("solved=$solved")
end # time_init_1

  if k == n
    if verb
      println("case k=n: solve directly")
    end
time_primal = @elapsed begin
    case = -2
    # xbarsort = x0sort - (s0 - r)/n
    lam = (s0 - r) / k
    @simd for i in 1:k
      @inbounds xbarsort[i] -= lam
    end
end # time_primal
    if hist
      return 0, (0,n), 0, (time_init_1, 0.0, time_primal), [(0,0)], case, lam
    else
      return 0, (0,n), 0, (time_init_1, 0.0, time_primal), case, lam
    end
  elseif k == 1
    if verb
      println("case k=1: solve directly")
    end
time_primal = @elapsed begin
    case = -1
    # xbarsort = min(x0sort, r)
    if x0prepop
      for i in 1:n
        @inbounds xbarsort[i] = min(x0sort[i], r)
        if x0sort[i] <= r
          break
        end
      end
    else
      @simd for i in 1:n
        @inbounds xbarsort[i] = min(x0sort[i], r)
      end
    end
end # time_primal
    if hist
      return 0, (0, findfirst(x0sort .< r)), 0, (time_init_1, 0.0, time_primal), [(0,0)], case, lam
    else
      return 0, (0, findfirst(x0sort .< r)), 0, (time_init_1, 0.0, time_primal), case, lam
    end
  end

time_init = @elapsed begin
  # preprocessing
  tol = zero(Tfr)
  if Tfr <: Union{Integer, Rational}
    tol = 0
  else
    # tol = eps(Tfr)*maximum(x0sort)
    tol = eps(Tfr)*x0sort[1]
  end
  case = 1
  n::Ti = length(x0sort)
  k0::Ti = k-1
  k1::Ti = k
  kk0::Ti = 0
  k1k0::Ti = 0
  the = zero(Tfr) # theta
  lam = zero(Tfr) # lambda
  theplam = zero(Tfr) # theta + lambda
  kkt2 = false
  kkt5 = false
  if verb
    kkt1 = false
    kkt3 = false
    kkt4 = false
  end
  nit::Ti = 0
  maxnit::Ti = length(x0sort)
  # sum1k0::Tfr = sum(view(x0sort, 1:k0))
  sum1k0::Tfr = s0 - x0sort[k]
  sumk0p1k1::Tfr = x0sort[k]
  sum1k0r = zero(Tfr)
  # sum1k0r_rho = zero(Tfr)
  # sumk0p1k1_rho = zero(Tfr)
  rho = zero(Tfr)
  solved = false
  if hist
    history = Vector{Tuple}(undef, n+1)
    history[1] = (k0, k1)
  end
end # time_init

  # iterate
time_pivot = @elapsed begin
  while true

    # counter
    nit += 1

    # rho, theta, and lambda
    kk0 = k - k0
    k1k0 = k1 - k0
    rho = k0 * k1k0 + kk0^2
    sum1k0r = sum1k0 - r
    the = (k0 * sumk0p1k1 - kk0 * sum1k0r) / rho #! OK
    # lam = (kk0 * sumk0p1k1 + k1k0 * (sum1k0 - r)) / rho #! OK
    #! OK
    if k0 > 0
      theplam = (k * the + sum1k0r) / k0
    else
      theplam = (k * sumk0p1k1 + (k1-k) * sum1k0r) / rho
    end

    # kkt conditions
    # kkt2 = (k0 == 0 ? true : x0sort[k0] > the + lam) #! OK
    # kkt2 = (k0 == 0 ? true : x0sort[k0] > theplam) #! OK
    # kkt5 = (k1 == n ? true : the > x0sort[k1+1]) #! OK
    kkt2 = (k0 == 0 ? true : x0sort[k0] > theplam-tol) #! OK
    kkt5 = (k1 == n ? true : the > x0sort[k1+1]-tol) #! OK
    
    # debug
    if verb
      kkt1 = lam > 0
      kkt3 = the + lam >= x0sort[k0+1]
      kkt4 = x0sort[k1] >= the
      println("----------------")
      println("k0=$k0, k1=$k1")
      # println("the=$the, lam=$lam, rho=$rho")
      println("the=$the, lam=$(theplam-the), rho=$rho")
      # println("sum1k0r_rho=$sum1k0r_rho, sumk0p1k1_rho=$sumk0p1k1_rho")
      println("sum1k0r_rho=$(sum1k0r/rho), sumk0p1k1_rho=$(sumk0p1k1/rho)")
      # println("kkt1=$kkt1, val=$(lam)")
      println("kkt1=$kkt1, val=$(theplam-the)")
      # println("kkt2=$kkt2, val=$(k0==0 ? Inf : x0sort[k0]-(the+lam))")
      if k0 > 0
        println("kkt2=$kkt2, val=$(k0==0 ? Inf : x0sort[k0]-(theplam))", ", $(x0sort[k0]) > $(theplam))")
      end
      # println("kkt3=$kkt3, val=$((the+lam) - x0sort[k0+1])")
      println("kkt3=$kkt3, val=$((theplam) - x0sort[k0+1])")
      println("kkt4=$kkt4, val=$(x0sort[k1] - the)")
      if k1 < n
        println("kkt5=$kkt5, val=$(the - (k1==n ? -Inf : x0sort[k1+1]))")
      end
    end

    # move k0
    if (kkt2 && kkt5) || (k0 == 0 && k1 == n)
      solved = true
      break
    elseif kkt2 && (k1 < n)
      k1 += 1
      sumk0p1k1 += x0sort[k1]
      # @assert(sumk0p1k1 == sum(x0sort[k0+1:k1]))
      if verb
        printstyled("k1 <- k1+1\n"; color=:green)
      end
    elseif !kkt2 && (k0 > 0)
      sum1k0 -= x0sort[k0]
      sumk0p1k1 += (k0 > 0 ? x0sort[k0] : 0)
      k0 -= 1
      # @assert(sumk0p1k1 == sum(x0sort[k0+1:k1]))
      if verb
        printstyled("k0 <- k0-1\n"; color=:blue)
      end
    end

    # history
    if hist
      history[nit+1] = (k0, k1)
    end
  end
end # time_pivot

  # stop
time_primal = @elapsed begin
  if solved
    lam = theplam - the
    if x0prepop
      @simd for i in 1:k0
        @inbounds xbarsort[i] -= lam
      end
    else
      @simd for i in 1:k0
        @inbounds xbarsort[i] = x0sort[i] - lam
      end
    end
    @simd for i in k0+1:k1
      @inbounds xbarsort[i] = the
    end
    if !x0prepop
      @simd for i in k1+1:n
        @inbounds xbarsort[i] = x0sort[i]
      end
    end
  end
end # time_primal

  if solved
    if !hist
      return 1, (k0, k1), nit, (time_init, time_pivot, time_primal), case, lam
    else
      return 1, (k0, k1), nit, (time_init, time_pivot, time_primal), history, case, lam
    end
  elseif nit >= maxnit
    if !hist
      if verb; printstyled("error\n"; color=:red); end
      return 1, (-2, -2), nit, (time_init, time_pivot, time_primal), case, lam
    else
      if verb; printstyled("error\n"; color=:red); end
      return 1, (-2, -2), nit, (time_init, time_pivot, time_primal), history, case, lam
    end
  end
end

#=
"""alternate numerical approaches for snake; alt and rho are wrong??"""
function project_maxksum_snakealt!(
  xbarsort::AbstractVector{Tfr}, x0sort::AbstractVector{Tfr}, r::Tfr, k::Ti,
  active::Bool, x0prepop::Bool=false, verb::Bool=false, hist::Bool=false,
) where {Tfr<:Union{AbstractFloat,Rational},Ti<:Integer}
  time_init = @elapsed begin
  
  # preprocessing
  tol = zero(Tfr)
  if Tfr <: Union{Integer, Rational}
    tol = 0
  else
    tol = eps(Tfr)*maximum(x0sort)
  end
  n::Ti = length(x0sort)
  k0::Ti = k-1
  k1::Ti = k
  kk0::Ti = 0
  k1k0::Ti = 0
  the = zero(Tfr)
  lam = zero(Tfr)
  kkt2 = false
  kkt5 = false
  if verb
    kkt1 = false
    kkt3 = false
    kkt4 = false
  end
  nit::Ti = 0
  maxnit::Ti = length(x0sort)
  sum1k0::Tfr = sum(x0sort[1:k0])
  sumk0p1k1::Tfr = sum(x0sort[k0+1:k1])
  sum1k1::Tfr = sum(x0sort[1:k1])
  rho = zero(Tfr)
  solved = false
  if hist
    history = Vector{Tuple}(undef, n+1)
    history[1] = (k0, k1)
  end
  end # time_init

  if !active
    if hist
      return 0, (-1, -1), nit, (time_init, 0.0, 0.0), hist
    else
      return 0, (-1, -1), nit, (time_init, 0.0, 0.0)
    end
  end

  # iterate
  time_pivot = @elapsed begin
  while true

    # counter
    nit += 1

    # rho, theta, and lambda
    #=
    #! "k/rho" idea
    rho = k0 * (k1-k0) + (k-k0)^2
    kk0_rho = (k - k0) / rho
    k1k0_rho = (k1 - k0) / rho
    the = (k0/rho) * sumk0p1k1 - kk0_rho * (sum1k0 - r)
    lam = kk0_rho * sumk0p1k1 + k1k0_rho * (sum1k0 - r)
    =#
    #! "s1k1, s1k0" idea
    rho = k0 * (k1-k0) + (k-k0)^2
    the = (k0/rho) * sum1k1 - (k/rho) * (sum1k0 - r)
    lam = (k/rho) * sumk0p1k1 + (k1/rho) * sum1k0 - (k0/rho) * sum1k1 + ((k0-k1)/rho) * r
    println("the=$the")
    println("lam=$lam\n")
    # the = (k0 * sum1k1 - k * (sum1k0 - r)) / rho
    the = (k0/rho) * sumk0p1k1 - (k-k0)/rho * (sum1k0 - r)
    lam = (k * sumk0p1k1 + k1 * sum1k0 - k0 * sum1k1 + (k0-k1) * r) / rho
    println("the=$the")
    println("lam=$lam\n")
    
    # kkt conditions
    # kkt2 = (k0 == 0 ? true : x0sort[k0] > the + lam + tol) #! dont use this; when are strict inequalities here an issue...
    # kkt5 = (k1 == n ? true : the > x0sort[k1+1] + tol) #! dont use this; when are strict inequalities here an issue...
    kkt2 = (k0 == 0 ? true : x0sort[k0] > the + lam)
    kkt5 = (k1 == n ? true : the > x0sort[k1+1])
    
    # debug
    if verb
      kkt1 = lam > 0
      kkt3 = the + lam >= x0sort[k0+1]
      kkt4 = x0sort[k1] >= the
      println("----------------")
      println("k0=$k0, k1=$k1")
      println("the=$the, lam=$lam, rho=$rho")
      println("kkt1=$kkt1, val=$(lam)")
      println("kkt2=$kkt2, val=$(k0==0 ? Inf : x0sort[k0]-(the+lam))")
      println("kkt3=$kkt3, val=$((the+lam) - x0sort[k0+1])")
      println("kkt4=$kkt4, val=$(x0sort[k1] - the)")
      println("kkt5=$kkt5, val=$(the - (k1==n ? -Inf : x0sort[k1+1]))")
    end

    # move k0
    if (kkt2 && kkt5) || (k0 == 0 && k1 == n)
      solved = true
      break
    elseif kkt2 && (k1 < n)
      k1 += 1
      sumk0p1k1 += x0sort[k1]
      sum1k1 += x0sort[k1]
      # @assert(sumk0p1k1 == sum(x0sort[k0+1:k1]))
      if verb
        printstyled("k1 <- k1+1\n"; color=:green)
      end
    elseif !kkt2 && (k0 > 0)
      sum1k0 -= x0sort[k0]
      sumk0p1k1 += (k0 > 0 ? x0sort[k0] : 0)
      k0 -= 1
      # @assert(sumk0p1k1 == sum(x0sort[k0+1:k1]))
      if verb
        printstyled("k0 <- k0-1\n"; color=:blue)
      end
    end
    if nit >= 2
      break
    end

    # history
    if hist
      history[nit+1] = (k0, k1)
    end
  end
  end # time_pivot

  # stop
  time_primal = @elapsed begin
  if solved
    if x0prepop
      @simd for i in 1:k0
        @inbounds xbarsort[i] = x0sort[i] - lam
      end
    else
      @simd for i in 1:k0
        @inbounds xbarsort[i] -= lam
      end
    end
    @simd for i in k0+1:k1
      @inbounds xbarsort[i] = the
    end
    if !x0prepop
      @simd for i in k1+1:n
        @inbounds xbarsort[i] = x0sort[i]
      end
    end
  end
  end # time_primal

  if solved
    if !hist
      return 1, (k0, k1), nit, (time_init, time_pivot, time_primal)
    else
      return 1, (k0, k1), nit, (time_init, time_pivot, time_primal), history
    end
  elseif nit >= maxnit
    if !hist
      if verb; printstyled("error\n"; color=:red); end
      return 1, (-2, -2), nit, (time_init, time_pivot, time_primal)
    else
      if verb; printstyled("error\n"; color=:red); end
      return 1, (-2, -2), nit, (time_init, time_pivot, time_primal), history
    end
  end
end

function project_maxksum_snakerho!(
  xbarsort::AbstractVector{Tfr}, x0sort::AbstractVector{Tfr}, r::Tfr, k::Ti,
  active::Bool, x0prepop::Bool=false, verb::Bool=false, hist::Bool=false,
) where {Tfr<:Union{AbstractFloat,Rational},Ti<:Integer}
  time_init = @elapsed begin
  
  # preprocessing
  tol = zero(Tfr)
  if Tfr <: Union{Integer, Rational}
    tol = 0
  else
    tol = eps(Tfr)*maximum(x0sort)
  end
  n::Ti = length(x0sort)
  k0::Ti = k-1
  k1::Ti = k
  kk0::Ti = 0
  k1k0::Ti = 0
  eta_the = zero(Tfr) # numerator of theta
  eta_lam = zero(Tfr) # numerator of lambda
  eta_theplam = zero(Tfr) # numerator of (theta+lambda)
  the = zero(Tfr)
  lam = zero(Tfr)
  kkt2 = false
  kkt5 = false
  if verb
    kkt1 = false
    kkt3 = false
    kkt4 = false
  end
  nit::Ti = 0
  maxnit::Ti = length(x0sort)
  sum1k0::Tfr = sum(x0sort[1:k0])
  sumk0p1k1::Tfr = sum(x0sort[k0+1:k1])
  sum1k0r = zero(Tfr)
  rho = zero(Tfr)
  solved = false
  if hist
    history = Vector{Tuple}(undef, n+1)
    history[1] = (k0, k1)
  end
  end # time_init

  if !active
    if hist
      return 0, (-1, -1), nit, (time_init, 0.0, 0.0), hist
    else
      return 0, (-1, -1), nit, (time_init, 0.0, 0.0)
    end
  end

  # iterate
  time_pivot = @elapsed begin
  while true

    # counter
    nit += 1

    # rho, theta, and lambda
    kk0 = k - k0
    k1k0 = k1 - k0
    rho = k0 * k1k0 + kk0^2
    sum1k0r = sum1k0 - r
    eta_the = k0 * sumk0p1k1 - kk0 * sum1k0r
    eta_lam = kk0 * sumk0p1k1 + k1k0 * sum1k0r
    eta_theplam = k * sumk0p1k1 + (k1-k) * sum1k0r
    if nit == 1
      eta_lam = 0.0007575390960099071 * rho
    end
    println("eta_the = k0 * sumk0p1k1 - kk0 * sum1k0r   = $(k0 * sumk0p1k1) - $(kk0 * sum1k0r)   = $eta_the")
    println("eta_lam = kk0 * sumk0p1k1 + k1k0 * sum1k0r = $(kk0 * sumk0p1k1) + $(k1k0 * sum1k0r) = $eta_lam")
    println("eta_tpl = k * sumk0p1k1 + k1k0 * sum1k0r   = $(k * sumk0p1k1) + $(k1k0 * sum1k0r)   = $eta_theplam")

    -9978.361931422378 + 7.567815569138972
    -9962.228222795024
    # kkt conditions
    # kkt2 = (k0 == 0 ? true : x0sort[k0] > the + lam + tol) #! dont use this; when are strict inequalities here an issue...
    # kkt5 = (k1 == n ? true : the > x0sort[k1+1] + tol) #! dont use this; when are strict inequalities here an issue...
    println("kkt2: rho * x0sort[k0] > eta_the + eta_lam == $(rho * x0sort[k0]) > $(eta_the + eta_lam)")
    println("kkt2: rho * x0sort[k0] > eta_theplam       == $(rho * x0sort[k0]) > $(eta_theplam)")
    println("kkt5: eta_the > rho * x0sort[k1+1] == $eta_the > $(rho * x0sort[k1+1])")
    # kkt2 = (k0 == 0 ? true : rho * x0sort[k0] > eta_the + eta_lam)
    kkt2 = (k0 == 0 ? true : rho * x0sort[k0] > eta_theplam)
    kkt5 = (k1 == n ? true : eta_the > rho * x0sort[k1+1])
    
    # debug
    if verb
      kkt1 = eta_lam > 0
      kkt3 = eta_the + eta_lam >= rho*x0sort[k0+1]
      kkt4 = rho*x0sort[k1] >= eta_the
      println("----------------")
      println("k0=$k0, k1=$k1")
      println("eta_the=$eta_the, eta_lam=$eta_lam, rho=$rho")
      println("sum1k0r=$sum1k0r, sumk0p1k1=$sumk0p1k1")
      println("kkt1=$kkt1, val=$(eta_lam)")
      println("kkt2=$kkt2, val=$(k0==0 ? Inf : rho*x0sort[k0]-(eta_the+eta_lam))")
      println("kkt3=$kkt3, val=$((eta_the+eta_lam) - rho*x0sort[k0+1])")
      println("kkt4=$kkt4, val=$(rho*x0sort[k1] - eta_the)")
      println("kkt5=$kkt5, val=$(eta_the - (k1==n ? -Inf : rho*x0sort[k1+1]))")
    end

    # move k0
    if (kkt2 && kkt5) || (k0 == 0 && k1 == n)
      solved = true
      break
    elseif kkt2 && (k1 < n)
      k1 += 1
      sumk0p1k1 += x0sort[k1]
      # @assert(sumk0p1k1 == sum(x0sort[k0+1:k1]))
      if verb
        printstyled("k1 <- k1+1\n"; color=:green)
      end
    elseif !kkt2 && (k0 > 0)
      sum1k0 -= x0sort[k0]
      sumk0p1k1 += (k0 > 0 ? x0sort[k0] : 0)
      k0 -= 1
      # @assert(sumk0p1k1 == sum(x0sort[k0+1:k1]))
      if verb
        printstyled("k0 <- k0-1\n"; color=:blue)
      end
    end

    # history
    if hist
      history[nit+1] = (k0, k1)
    end
  end
  end # time_pivot

  # stop
  time_primal = @elapsed begin
  the = eta_the/rho
  lam = eta_lam/rho
  if solved
    if x0prepop
      @simd for i in 1:k0
        @inbounds xbarsort[i] = x0sort[i] - lam
      end
    else
      @simd for i in 1:k0
        @inbounds xbarsort[i] -= lam
      end
    end
    @simd for i in k0+1:k1
      @inbounds xbarsort[i] = the
    end
    if !x0prepop
      @simd for i in k1+1:n
        @inbounds xbarsort[i] = x0sort[i]
      end
    end
  end
  end # time_primal

  if solved
    if !hist
      return 1, (k0, k1), nit, (time_init, time_pivot, time_primal)
    else
      return 1, (k0, k1), nit, (time_init, time_pivot, time_primal), history
    end
  elseif nit >= maxnit
    if !hist
      if verb; printstyled("error\n"; color=:red); end
      return 1, (-2, -2), nit, (time_init, time_pivot, time_primal)
    else
      if verb; printstyled("error\n"; color=:red); end
      return 1, (-2, -2), nit, (time_init, time_pivot, time_primal), history
    end
  end
end
=#