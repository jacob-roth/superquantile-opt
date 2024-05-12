include("projection_lcp.jl")
include("projection_snake.jl")
include("projection_grid.jl")
include("projection_gurobi.jl")

"""
write solution of (`l`, `u`)-box projection of `x0` into `x`, where `l` <= `u` is NOT checked
"""
function project_box!(x::AbstractVector{Tf}, x0::AbstractVector{Tf}, l::AbstractVector{Tf}, u::AbstractVector{Tf}) where {Tf<:AbstractFloat}
  x .= min.(max.(x0, l), u)
  nothing
end
function freeidx_box!(j::AbstractVector{Tf}, x::AbstractVector{Tf}, l::AbstractVector{Tf}, u::AbstractVector{Tf}, tol::Tf=1e-10) where {Tf<:AbstractFloat}
  j .= (x .<= u .+ tol) .&& (x .>= l .- tol)
  nothing
end
# function considx_box!(j::AbstractVector{Tf}, x::AbstractVector{Tf}, l::AbstractVector{Tf}, u::AbstractVector{Tf}, tol::Tf=1e-10) where {Tf<:AbstractFloat}
#   j .= (x .>= u .+ tol) .|| (x .<= l .- tol)
#   nothing
# end
function considx_box!(j::AbstractVector{Tf}, x0::AbstractVector{Tf}, l::AbstractVector{Tf}, u::AbstractVector{Tf}, tol::Tf=1e-10) where {Tf<:AbstractFloat}
  # j .= (x0 .>= u .- tol) .|| (x0 .<= l .+ tol) # `tol` buffer
  j .= (x0 .>= u .+ tol) .|| (x0 .<= l .- tol) # allow `tol` infeasibility
  nothing
end
