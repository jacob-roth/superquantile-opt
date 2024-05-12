struct ProjectBWorkspace{
  Tf<:AbstractFloat,Vf<:AbstractMatrix{Tf},VVf<:AbstractVector{Vf},
  Ti<:Integer,Vi<:AbstractMatrix{Ti},VVi<:AbstractVector{Vi}
}
  ytilde::VVf # input
  ybar::VVf # projection
  yhat::VVf # ytilde - ybar
  sig::VVi # ordering
  ytmp::VVf # temporary
end

function project_B!(
  ytilde::VVf, sig::VVi, # already updated
  ybar::VVf, yhat::VVf, ytmp::VVf, # ProjBWorkspace
  k0bar::Vi, k1bar::Vi, # diagnostics
  r::Vf, k::Vi, active::Vector{Bool}, redo::Vector{Bool}, # parameters
) where {
  Tf<:AbstractFloat,Vf<:AbstractVector{Tf},VVf<:AbstractVector{Vf},
  Ti<:Integer,Vi<:AbstractVector{Ti},VVi<:AbstractVector{Vi}
}
  eta::Tf = 0.0
  sp::Ti = 0
  spl = 0
  k0l = 0
  k1l = 0
  for l in eachindex(ytilde)
    if !redo[l]
      continue
    end

    @views eta = norm(ytilde[l], Inf)
    @views ytmp[l] .= ytilde[l][sig[l]] ./ eta #! bottleneck
    spl, (k0l, k1l), _ = @views project_maxksum_snake!(ybar[l], ytmp[l], r[l] / eta, k[l], active[l], false) #! here ybar[l] is assumed ("incorrectly") to be sorted
    sp += spl
    @views ytmp[l] .= ybar[l] .* eta
    @views ybar[l][sig[l]] .= ytmp[l] #! sort ybar[l]
    @views yhat[l] .= ytilde[l] .- ybar[l];
    
    # k0bar, k1bar
    k0bar[l] = k0l
    k1bar[l] = k1l
  end
  return sp
end

struct ProjectXWorkspace{Tf<:AbstractFloat,Vf<:AbstractVector{Tf}}
  ztilde::Vf # input
  zbar::Vf # projection
  zhat::Vf # ztilde - zbar
end

function project_X!(
  ztilde::Vf, zbar::Vf, zhat::Vf, # ProjXWorkspace
  xlo::Vf, xhi::Vf # parameters
) where {Tf<:AbstractFloat,Vf<:AbstractVector{Tf}}
  project_box!(zbar, ztilde, xlo, xhi)
  zhat .= ztilde .- zbar
  nothing
end