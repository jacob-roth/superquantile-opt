function pobj!(
  FV::FunctionValue, x::Vf, C::NDf, c::Vf, Cx::Vf,
) where {
  Tf<:AbstractFloat,Vf<:AbstractVector{Tf},Df<:Diagonal{Tf},NDf<:Union{Nothing,Df}
}
  if !isnothing(C)
    Cx .= C * x
    FV.v = 0.5dot(Cx, x) + dot(c, x)
    FV.g .= Cx .+ c
  else
    FV.v = dot(c, x)
  end
  return FV.v
end

function dobj!(
  FV::FunctionValue, objtype::Ti,
  lambda::VVf, mux::Vf,
  negmux_neg::Vf, negmux_pos::Vf,
  Atlampmu::Vf, Btlam::Vector{Vf},
  b::Vector{Vf}, r::Vector{Tf}, k::Vector{Ti}, 
  xlo::Vf, xhi::Vf, xloidx::Vector{Bool}, xhiidx::Vector{Bool}, Xunconstr::Bool,
  Cinv::NDf, c::Vf,
  tol::Tf=1e-8,
) where {
  Tf<:AbstractFloat,Ti<:Integer,
  Vf<:AbstractVector{Tf},
  NDf<:Union{Nothing,Diagonal{Tf}},
  VVf<:AbstractVector{Vf}
}
  # if prod(any(Atildeinv[l]' * lambda[l] .>= tol) for l in eachindex(lambda))
  bdd = true
  for l in eachindex(lambda)
    @inbounds for i in eachindex(lambda[l])
      bdd *= lambda[l][i] <= tol # (-lambda) >= 0: (-lambda) >= -tol <==> lambda <= tol
      if bdd == false
        break
      end
    end
  end
  if !bdd
    FV.d = -Inf
  # elseif prod(any(xhi .< xlo)) #!handled elsewhere!
  #    return +Inf
  else
    if objtype == 2
      FV.d = -0.5((Atlampmu - c)' * Cinv * (Atlampmu - c)) - b'*lambda + sum((r[l]/k[l]) * sum(lambda[l]) for l in eachindex(lambda))
    else
      FV.d = -b'*lambda + sum((r[l]/k[l]) * sum(lambda[l]) for l in eachindex(lambda))
    end
    if !Xunconstr
      negmux_neg .= min.(-mux, 0.0)
      negmux_pos .= max.(-mux, 0.0)
      @views FV.d -= (xlo[xloidx]' * negmux_neg[xloidx] + xhi[xhiidx]' * negmux_pos[xhiidx])
      # FV.d -= (xlo' * min.(-mux, 0.0) + xhi' * max.(-mux, 0.0))
    end
  end
end
