function vphi!(
  phi::FV1, obj::FV2,
  x::Vf, xprev::Vf,
  yhat::Vector{Vf},
  uhat::Vf,
  sigma::Tf, tau::Tf, PI::ProblemInstance,
  Cx::Union{Nothing,Vf}, BXratio::Tf=1.0,
) where {FV1<:FunctionValue, FV2<:FunctionValue, Tf<:AbstractFloat, Vf<:AbstractVector{Tf}}
ctime = @elapsed begin
  # evaluate phi
  pobj!(obj, x, PI.C, PI.c, Cx)
  if !PI.Xunconstr
    phi.v = obj.v + sigma/2 * norm(yhat, 2)^2 + (sigma * BXratio)/2 * norm(uhat, 2)^2
  else
    phi.v = obj.v + sigma/2 * norm(yhat, 2)^2
  end
  
  # prox adjustment
  if tau > 0
     phi.v += tau/(2*sigma) * norm(x .- xprev, 2)^2
  end
end # ctime
  return ctime
end
function gphi!(
  phi::FV1, obj::FV2,
  x::Vf, xprev::Vf,
  yhat::Vector{Vf},
  uhat::Vf,
  sigma::Tf, tau::Tf, PI::ProblemInstance,
  Cx::Union{Nothing,Vf}, BXratio::Tf=1.0,
) where {FV1<:FunctionValue, FV2<:FunctionValue, Tf<:AbstractFloat, Vf<:AbstractVector{Tf}}
ctime_gphi = @elapsed begin  
  # evaluate ∇phi
  pobj!(obj, x, PI.C, PI.c, Cx)
  phi.g .= obj.g
  for l in 1:PI.L
    # phi.g .+= sigma .* (PI.A[l]' * (yhat[l]))
    # compute via BLAS: phi.g .+= sigma .* (PI.A[l]' * (yhat[l]))
    # gemv!(tA, alpha, A, x, beta, y) => y = alpha*A*x + beta*y
    # @views BLAS.gemv!('T', sigma, PI.A[l], yhat[l], true, phi.g)
    
    # @views mul!(phi.g, PI.A[l]', yhat[l], sigma, true)
    if length(sparse(yhat[l]).nzind)/length(yhat[l]) <= 0.02
      @views mul!(phi.g, PI.A[l]', sparse(yhat[l]), sigma, true) # O(k1bar*n)
    else
      @views mul!(phi.g, PI.A[l]', yhat[l], sigma, true) # O(m*n)
    end
  end
  if !PI.Xunconstr
    phi.g .+= (sigma * BXratio) .* (uhat)
  end
  
  # prox adjustment
  if tau > 0
    phi.g .+= (tau/sigma) .* (x .- xprev)
  end
  # println("||g|| = $(norm(phi.g))")
end # ctime_gphi
  return ctime_gphi
end
function prepareHphi!(
  x::Vf, xprev::Vf, d::Vf,
  ytilde::Vector{Vf}, ybar::Vector{Vf}, ytmp::Vector{Vf}, sig::Vector{Vi}, k0bar::Vi, k1bar::Vi,
  utilde::Vf, ubar::Vf,
  sigma::Tf, tau::Tf,
  phi::FV1, obj::FV2, PI::ProblemInstance,
  Ttilde::Mf, Ttilde_rowmax::Base.RefValue{Ti}, c::Vector{Vr}, a::Vector{Vf}, b::Vector{Vf}, active::Vector{Bool}, diagjX::Vf,
  debug::Bool=false, BXratio::Tf=1.0,
) where {
  Tf<:AbstractFloat,Vf<:AbstractVector{Tf},Mf<:AbstractMatrix{Tf},
  Ti<:Integer,Vi<:AbstractVector{Ti},
  Tr<:Rational{Ti},Vr<:AbstractVector{Tr},
  FV1<:FunctionValue,FV2<:FunctionValue
}
  #  
  # initialize
  #
  
  idx_lo = 0
  idx_hi = 0
  idx = Vector{UnitRange}(undef, PI.L)
  
  #
  # max-k-sum
  #

  for l in 1:PI.L
    if active[l]
      # if k1bar[l] == PI.k[l]; @warn("l=$l, k1bar == k"); end
      
      if k1bar[l] != PI.k[l]# Case 2
      # if true # ONLY k=10% L=10. obj=1 run with this!
        # a, b, c, and Ttilde
        idx_lo = idx_hi+1
        idx_hi = idx_lo+(k1bar[l]-k0bar[l])
        idx[l] = idx_lo:idx_hi
        @views set_c!(c[l], PI.k[l], k0bar[l], k1bar[l])
        #? idea: single for-loop over rows in Ttilde?
        #? idea: G[siginv,siginv] * v = G[siginv,:] * v[sig]
        @views Arows1 = PI.A[l][sig[l][1:k0bar[l]],:]
        @views Arows2 = PI.A[l][sig[l][k0bar[l]+1:k1bar[l]],:]
        @views a[l] .= sum(Arows1, dims=1)' #! bottleneck
        @views b[l] .= sum(Arows2, dims=1)' #! bottleneck
        @views Ttilde[idx_lo,:] .= ((sqrt(k0bar[l]) * c[l][1]) .* a[l] .+ (sqrt(k0bar[l]) * c[l][2]) .* b[l])
        # @views Ttilde[idx_lo+1:idx_hi,:] .= Arows2 .+ c[l][2] .* a[l]' .- c[l][3] .* b[l]'  #! bottleneck
        @views Ttilde[idx_lo+1:idx_hi,:] .= Arows2 .+ c[l][2] .* a[l]' .- c[l][3] .* b[l]'  #! bottleneck
      else# Case 3
        println("CASE 3")
        idx_lo = idx_hi+1
        idx_hi = idx_lo
        idx[l] = idx_lo:idx_hi
        # println("IDX=$(idx[l])")
        # println("dim(Ttilde)=$(size(Ttilde[idx[l],:]))")
        @views Arows_alphabeta = PI.A[l][sig[l][1:k1bar[l]],:]
        # println("dim fill = $(size(sum(Arows_alphabeta, dims=1)'/sqrt(k1bar[l])))")
        @views Ttilde[idx[l],:] .= sum(Arows_alphabeta, dims=1)/sqrt(k1bar[l])
      end
    end
  end

  if debug
    println("========")
    println("idx hi: ", idx_hi)
    println("========")
  end
  
  if idx_hi == 0
    Ttilde .= 0.0
  elseif idx_hi < Ttilde_rowmax[]
    #!!!!!!! idx_hi+1:end !!!!!!!
    #! bottleneck
    @views Ttilde[idx_hi+1:Ttilde_rowmax[],:] .= 0.0
  end
  Ttilde_rowmax[] = max(Ttilde_rowmax[], idx_hi)
  
  #
  # box constraint
  #
  
  if !PI.Xunconstr
    considx_box!(diagjX, utilde, PI.xlo, PI.xhi, eps(eltype(utilde)))
  end
  return idx_lo, idx_hi
end

function Hphi!(
  x::Vf, xprev::Vf, d::Vf,
  ytilde::Vector{Vf}, ybar::Vector{Vf}, ytmp::Vector{Vf}, sig::Vector{Vi}, k0bar::Vi, k1bar::Vi,
  utilde::Vf, ubar::Vf,
  sigma::Tf, tau::Tf,
  phi::FV1, obj::FV2, PI::ProblemInstance,
  Cx::Union{Nothing,Vf},
  Ttilde::Mf, Ttilde_rowmax::Base.RefValue{Ti}, c::Vector{Vr}, a::Vector{Vf}, b::Vector{Vf}, active::Vector{Bool}, diagjX::Vf,
  debug::Bool=false, BXratio::Tf=1.0,
) where {
  Tf<:AbstractFloat,Vf<:AbstractVector{Tf},Mf<:AbstractMatrix{Tf},
  Ti<:Integer,Vi<:AbstractVector{Ti},
  Tr<:Rational{Ti},Vr<:AbstractVector{Tr},
  FV1<:FunctionValue,FV2<:FunctionValue
}
ctime = @elapsed begin

  pobj!(obj, x, PI.C, PI.c, Cx)
  
  if sum(active) == 0
    phi.H .= obj.H
    if tau > 0
      t_s = tau/sigma
      for i in 1:PI.n
        phi.H[i,i] += t_s
      end
    end
  else
    idx_lo, idx_hi = prepareHphi!(
      x, xprev, d,
      ytilde, ybar, ytmp, sig, k0bar, k1bar,
      utilde, ubar,
      sigma, tau,
      phi, obj, PI,
      Ttilde, Ttilde_rowmax, c, a, b, active, diagjX,
      debug, BXratio
    )
    # phi.H .= obj.H .+ sigma .* Ttilde[1:idx_hi,:]' * Ttilde[1:idx_hi,:] #! here
    phi.H .= obj.H
    @views mul!(phi.H, Ttilde[1:idx_hi,:]', Ttilde[1:idx_hi,:], sigma, true)
    # @time phi.H .= obj.H .+ sigma .* Ttilde[1:idx_hi,:]' * Ttilde[1:idx_hi,:] #! here
    # println("^^ T'*T")
    # prox adjustment
    if tau > 0
      t_s = tau/sigma
      for i in 1:PI.n
        phi.H[i,i] += t_s
      end
    end
  end
  
  # box constraint
  if !PI.Xunconstr
    phi.H .+= (sigma * BXratio) .* Diagonal(diagjX)
  end
end # ctime
  return ctime
end
function vgphi!(
  phi::FV1, obj::FV2,
  x::Vf, xprev::Vf,
  yhat::Vector{Vf},
  uhat::Vf,
  sigma::Tf, tau::Tf, PI::ProblemInstance,
  Cx::Union{Nothing,Vf}, BXratio::Tf=1.0,
) where {FV1<:FunctionValue, FV2<:FunctionValue, Tf<:AbstractFloat, Vf<:AbstractVector{Tf}}
  ctime_vphi = @elapsed begin
  # evaluate v phi
  pobj!(obj, x, PI.C, PI.c, Cx)
  phi.v = obj.v
  if !PI.Xunconstr
    phi.v += sigma/2 * norm(yhat, 2)^2 + (sigma * BXratio)/2 * norm(uhat, 2)^2
  else
    phi.v += sigma/2 * norm(yhat, 2)^2
  end
  end # ctime_vphi
  
  # evaluate g phi
  ctime_gphi = @elapsed begin
  phi.g .= obj.g
  for l in eachindex(PI.m)
    # phi.g .+= sigma .* (PI.A[l]' * (yhat[l]))
    # compute via BLAS: phi.g .+= sigma .* (PI.A[l]' * (yhat[l]))
    # gemv!(tA, alpha, A, x, beta, y) => y = alpha*A*x + beta*y
    # @views BLAS.gemv!('T', sigma, PI.A[l], yhat[l], true, phi.g)
    
    # prev
    # @views mul!(phi.g, PI.A[l]', yhat[l], sigma, true) # O(m*n)

    # new
    if length(sparse(yhat[l]).nzind)/length(yhat[l]) <= 0.02
      @views mul!(phi.g, PI.A[l]', sparse(yhat[l]), sigma, true) # O(k1bar*n)
    else
      @views mul!(phi.g, PI.A[l]', yhat[l], sigma, true) # O(m*n)
    end
  end
  if !PI.Xunconstr
    phi.g .+= (sigma * BXratio) .* (uhat)
  end
  # prox adjustment
  if tau > 0
    phi.v += tau/(2*sigma) * norm(x .- xprev, 2)^2
    phi.g .+= (tau/sigma) .* (x .- xprev)
  end
  # println("||g|| = $(norm(phi.g))")
  end # ctime_vphi
  return ctime_vphi, ctime_gphi
end
