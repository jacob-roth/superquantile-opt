function calc_ADMM_line(
  x::Vf, xprev::Vf, y::Vector{Vf}, z::Vf, lambda::Vector{Vf}, mu::Vf, negmux_neg::Vf, negmux_pos::Vf,
  Ax::Vector{Vf}, Btlam::Vector{Vf}, Atlampmu::Vf,
  obj::FunctionValue, sigma::Tf, tau::Tf,
  PI::ProblemInstance
) where {Tf <: AbstractFloat, Vf <: AbstractVector{Tf}}
  # calc line items
  pinfeas = maximum([[norm(y[l] .- PI.A[l]*x .- PI.b[l],Inf) for l in 1:PI.L]; norm(x .- z, Inf)])
  dinfeas = norm(obj.g .- (PI.A' * lambda + mu), Inf)
  pobj!(obj, x, PI.C, PI.c, zeros(PI.n))
  dobj!(obj, PI.objtype,
    lambda, mu,
    negmux_neg, negmux_pos,
    Atlampmu, Btlam, PI.b, PI.r, PI.k,
    PI.xlo, PI.xhi, PI.xloidx, PI.xhiidx, PI.Xunconstr,
    PI.Cinv, PI.c
  )
  pobj = obj.v
  dobj = obj.d
  if !isinf(dobj)
    relgap = (pobj - dobj) / (1 + abs(pobj) + abs(dobj))
  else
    dobj = -9.999e99
    relgap = +9.999e99
  end
  return (pinfeas, dinfeas, relgap, pobj, dobj, sigma)
end

function ADMM_line(
  iter::Integer,
  pinfeas::Float64, dinfeas::Float64, relgap::Float64,
  pobj::Float64, dobj::Float64,
  sigma::Float64, tau::Float64,
  id::Integer,
  etime::Float64,
  header::Bool=false,
)
  line = "";
  if header
    header = "";
    header *= "ALMit   "
    header *= "pinfeas      "
    header *= "dinfeas      "
    header *= "relgap       "
    header *= "pobj         "
    header *= "dobj        "
    header *= "sigma    tau        id   "
    # header *= "lsit        "
    header *= "etime\n"
    line *= header
  end
  line *= @sprintf("%5.0d", iter);
  line *= @sprintf("   %+9.3e", pinfeas);
  line *= @sprintf("   %+9.3e", dinfeas);
  line *= @sprintf("   %+9.3e", relgap);
  line *= @sprintf("   %+9.3e", pobj);
  line *= @sprintf("   %+9.3e", dobj);
  if id == 1
    line *= @sprintf("  %6.1e  %6.1e %6s", sigma, tau, "avg");
  elseif id == 2
    line *= @sprintf("  %6.1e  %6.1e %6s", sigma, tau, "itr");
  end
  # line *= @sprintf("   %9.0d", ls);
  line *= @sprintf("  %9.2e\n", etime);
  return line
end

function calc_ALM_line(
  x::Vf, xprev::Vf, y::Vector{Vf}, u::Vf,
  lambda::Vector{Vf}, mux::Vf, negmux_neg::Vf, negmux_pos::Vf,
  Ax::Vector{Vf}, Btlam::Vector{Vf},
  Atlampmu::Vf,
  phi::FunctionValue, obj::FunctionValue, sigma::Tf, tau::Tf, PI::ProblemInstance,
  Cx::Vf, Ry::Vector{Vf}, Ru::Vf, Rd::Vf, sig::Vector{Vi}
) where {Ti <: Integer, Tf <: AbstractFloat, Vf <: AbstractVector{Tf}, Vi <: AbstractVector{Ti}}
  # pinfeas = maximum([[norm(Ry[l], Inf) for l in 1:PI.L]; norm(Ru, Inf)])
  # dinfeas = norm(obj.g .- Atlampmu, Inf)
  for l in eachindex(PI.m)
    Ry[l] .= max.(Ax[l] .+ PI.b[l] .- y[l], 0.0) # y >= Ax + b
    # Ry[l] .= (Ax[l] .+ PI.b[l] .- y[l]) ./ (1+norm(PI.b[l])) # y == Ax + b
  end
  Ru .= u .- x
  Rd .= (obj.g .- Atlampmu) ./ (1+norm(obj.g, Inf))
  y_Axpb = maximum(norm.(Ry, Inf)./(1 .+ norm.(PI.b, Inf)))
  # y_Axpb = maximum([norm(Ry[l][sig[l][1:PI.k[l]]], 1) for l in eachindex(PI.m)])
  u_x = norm(Ru, Inf) / (1+norm(u))
  pinfeas = max(
    y_Axpb, u_x,
    maximum([max((sum((Ax[l] .+ PI.b[l])[sig[l][1:PI.k[l]]]) - PI.r[l])/PI.k[l], 0.0) for l in 1:PI.L])
  )
  dinfeas = norm(Rd, Inf)
  pobj!(obj, x, PI.C, PI.c, Cx)
  dobj!(obj, PI.objtype,
    lambda, mux,
    negmux_neg, negmux_pos,
    Atlampmu, Btlam, PI.b, PI.r, PI.k,
    PI.xlo, PI.xhi, PI.xloidx, PI.xhiidx, PI.Xunconstr,
    PI.Cinv, PI.c
  )
  pobj = obj.v
  dobj = obj.d
  if !isinf(dobj)
    # relgap = (pobj - dobj) / (1 + abs(pobj) + abs(dobj))
    relgap = (pobj - dobj) / (1 + abs(pobj))
  else
    dobj = -9.999e99
    relgap = +9.999e99
  end
  return (pinfeas, y_Axpb, u_x, dinfeas, relgap, pobj, dobj, sigma)
end

function ALM_line(
  iter::Real,
  pinfeas::Float64, pinfeas_goal::Float64, dinfeas::Float64, dinfeas_goal::Float64, relgap::Float64,
  pobj::Float64, dobj::Float64,
  sigma::Float64, tau::Float64,
  ls::Integer, nit::Integer, phase::Integer, tolsub::Float64, suberr::Float64,
  etime::Float64,
  ctime::Float64,
  header::Bool=false,
)
  line = "";
  if header
    header = "";
    header *= "ALMit   "
    header *= "pinfeas      "
    header *= "pinfeas_goal "
    header *= "dinfeas      "
    header *= "dinfeas_goal "
    header *= "relgap       "
    header *= "pobj         "
    header *= "dobj        "
    # header *= "sigma    tau       lsit    nit   "
    header *= "sigma    tau        nit   lsit  phase  tolsub   suberr    "
    # header *= "lsit        "
    header *= "etime      "
    header *= "ctime\n"
    line *= header
  end
  line *= @sprintf("%5.1f", iter);
  line *= @sprintf("   %+9.3e", pinfeas);
  line *= @sprintf("   %+9.3e", pinfeas_goal);
  line *= @sprintf("   %+9.3e", dinfeas);
  line *= @sprintf("   %+9.3e", dinfeas_goal);
  line *= @sprintf("   %+9.3e", relgap);
  line *= @sprintf("   %+9.3e", pobj);
  line *= @sprintf("   %+9.3e", dobj);
  # line *= @sprintf("  %6.1e  %6.1e %6.0d %6.0d", sigma, tau, ls, nit);
  line *= @sprintf("  %6.1e  %6.1e %6.0d %6.0d %6.0d  %6.1e  %6.1e", sigma, tau, nit, ls, phase, tolsub, suberr);
  # line *= @sprintf("   %9.0d", ls);
  line *= @sprintf("  %9.2e", etime);
  line *= @sprintf("  %9.2e\n", ctime);
  return line
end

function calc_SUB_line(
  x::Vf, xprev::Vf, y::Vector{Vf}, u::Vf, 
  lambda::Vector{Vf}, mux::Vf, negmux_neg::Vf, negmux_pos::Vf,
  Ax::Vector{Vf}, Btlam::Vector{Vf},
  Atlampmu::Vf,
  phi::FunctionValue, obj::FunctionValue, sigma::Tf, tau::Tf, PI::ProblemInstance,
  Cx::Vf, Ry::Vector{Vf}, Ru::Vf, Rd::Vf, sig::Vector{Vi}
) where {Ti <: Integer, Tf <: AbstractFloat, Vf <: AbstractVector{Tf}, Vi <: AbstractVector{Ti}}
  # pinfeas = maximum([[norm(y[l] .- Ax[l] .- PI.b[l],Inf) for l in 1:PI.L]; norm(x .- u, Inf)])
  # dinfeas = norm(obj.g .- Atlampmu, Inf)
  for l in eachindex(PI.m)
    # Ry[l] .= max.(Ax[l] .+ PI.b[l] .- y[l], 0.0) # y >= Ax + b
    Ry[l] .= (Ax[l] .+ PI.b[l] .- y[l]) # y == Ax + b
  end
  Ru .= u .- x
  Rd .= (obj.g .- Atlampmu) ./ (1+norm(obj.g, Inf))
  y_Axpb = maximum(norm.(Ry, Inf)./(1 .+ norm.(PI.b, Inf)))
  # y_Axpb = maximum([norm(Ry[l][sig[l][1:PI.k[l]]], 1) for l in eachindex(PI.m)])
  u_x = norm(Ru, Inf) / (1+norm(u))
  # pinfeas = max(y_Axpb, u_x)
  pinfeas = max(
    y_Axpb, u_x,
    maximum([max((sum((Ax[l] .+ PI.b[l])[sig[l][1:PI.k[l]]]) - PI.r[l])/PI.k[l], 0.0) for l in 1:PI.L])
  )
  dinfeas = norm(Rd, Inf)
  
  pobj!(obj, x, PI.C, PI.c, Cx)
  #! already computed dobj in outer
  # dobj!(obj, PI.objtype,
  #   lambda, mux,
  #   negmux_neg, negmux_pos,
  #   Atlampmu, Btlam, PI.b, PI.r, PI.k,
  #   PI.xlo, PI.xhi, PI.xloidx, PI.xhiidx, PI.Xunconstr,
  #   PI.Cinv, PI.c
  # )
  pobj = obj.v
  dobj = obj.d
  if !isinf(dobj)
    # relgap = (pobj - dobj) / (1 + abs(pobj) + abs(dobj))
    relgap = (pobj - dobj) / (1 + abs(pobj))
  else
    dobj = -9.999e99
    relgap = +9.999e99
  end

  # println("---------- diagnostics ----------")
  # println("  dobj = $dobj")
  # println("  pobj = $pobj")
  # println("  lam  = $lambda")
  # println("  b  = $(PI.b)")
  # println("  b*lam = $(lambda' * PI.b)")
  return (pinfeas, y_Axpb, u_x, dinfeas, pobj, relgap, phi.v, norm(phi.g))
end

function SUB_line(
  iter::Integer,
  pinfeas::Float64, y_Axpb::Float64, u_x::Float64, dinfeas::Float64,
  pobj::Float64, relgap::Float64,
  phix::Float64, normgphix::Float64,
  linsys::Int64, absbeta::Int64,
  stepl::Float64, lsit::Int64,
  active::Vector{Bool},
  etime::Float64,
  top::Bool, bot::Bool,
)
  factive = findall(active)
  active_print = ["", "", "", "  "]
  for i in eachindex(factive)
    active_print[min(i, 3)] = @sprintf("%d",factive[i])
  end
  if length(factive) >=4
    if factive[end] >= 100
      active_print[4] = "*"
    else
      active_print[4] = "* "
    end
  end
  nactive = length(factive)
  cols = ["pinfeas", "y==Ax+b", "u=x", "dinfeas", "relgap", "pobj", "φ(x)", "||∇φ(x)||", "linsysdim", "|β|",   "stepl", "lsit", "active", "nactive", "etime"]
  # cols = ["pinfeas", "y>=Ax+b", "u=x", "dinfeas", "relgap", "pobj", "φ(x)", "||∇φ(x)||", "linsysdim", "|β|",   "stepl", "lsit", "active", "nactive", "etime"]
  vals = [pinfeas,    y_Axpb,   u_x,   dinfeas,   relgap,   pobj,   phix,   normgphix,   linsys,      absbeta, stepl,   lsit, active_print,   nactive,   etime]
  @assert(length(cols) == length(vals))
  line = "";
  if top
    top = " ┌───────"
    for i in eachindex(cols)
      top *= "┬" * "─"^12
    end
    top *= "┐\n"
    
    header = ""
    header *= " │ SSNit │"
    for i in eachindex(cols)
      header *= @sprintf(" %-10s │", cols[i])
    end
    header *= "\n"
    line *= top * header
  end
  line *= iter == 0 ? @sprintf(" │ %5s ", "--") : @sprintf(" │ %5.0d ", iter);
  for i in eachindex(vals)
    if cols[i] == "active"
      if nactive == 0
        line *= @sprintf("│ %+10s ", "--")
      else
        line *= @sprintf("│ %2s %2s %2s %s", vals[i]...) #! only works for L<=99
      end
    elseif (cols[i] == "linsysdim" || cols[i] == "stepl" || cols[i] == "lsit" || cols[i] == "etime") && iter == 0
      line *= @sprintf("│ %+10s ", "--")
    else
      line *= @sprintf("│ %+9.3e ", vals[i])
    end
  end
  line *= "│\n"
  if bot
    bot = " ├───────"
    for i in eachindex(cols)
      bot *= "┴" * "─"^12
    end
    bot *= "┤\n"
    line *= bot
  end  
  return line
end

function SUB_line_term(dinfeas::Float64, normgphi::Float64, tolsub::Float64, good::Bool)
  line = ""
  if good
    line *= " └────────── GOOD subproblem termination: dinfeas = " * @sprintf("%+0.3e", dinfeas)
  else
    line *= " └─────────── BAD subproblem termination: dinfeas = " * @sprintf("%+0.3e", dinfeas)
  end
  line *= ", ||∇φ(x)|| = " * @sprintf("%+0.3e", normgphi)
  line *= ", tolsub = " * @sprintf("%+0.3e", tolsub) * " ────────────────────────────────────────────────────────────────────────────────────────────────┘\n"
  return line
end
  
#=
iter = 2345
pinfeas = 1.234e-3
dinfeas = 4.321e-3
relgap = 5.678e-4
pobj = 3.456
dobj = 3.456
sigma = 1.234e-2
ls = 1234
etime = 0.0123
a = ALM_line(iter, pinfeas, dinfeas, relgap, pobj, dobj, sigma, ls, etime)
s1 = SUB_line(1, pinfeas, dinfeas, relgap, pobj, dobj, sigma, 1, etime, true, false)
s2 = SUB_line(123, pinfeas, dinfeas, relgap, pobj, dobj, sigma, 12, etime, false, false)
s3 = SUB_line(12345, pinfeas, dinfeas, relgap, pobj, dobj, sigma, 123, etime, false, true)
t = SUB_line_term(1.03, 1.03, 1.03)
println(a * s1 * s2 * s3 * t)

iter = 23
pinfeas = -1.2e-3
dinfeas = 4.1e-3
relgap = -5.8e-4
pobj = 3.6
dobj = 3.6
etime = 0.01
sigma = 1.4e-2
ls = 1

iter = 23
pinfeas = -1.2e-3
dinfeas = -4.1e-3
relgap = -5.8e-4
pobj = -3.6
dobj = -3.6
etime = 0.01
sigma = 1.4e-2
ls = 1

SUB_line(iter, pinfeas, dinfeas, relgap, pobj, dobj, sigma, ls, etime)
********************************************************************************************************

ALMit   pinfeas      dinfeas      relgap       pobj         dobj         sigma        lsit   etime
   23   -1.200e-03   -4.100e-03   -5.800e-04   -3.600e+00   -3.600e+00   +1.400e-02      1   1.00e-02
ALMit   pinfeas      dinfeas      relgap       pobj         dobj         sigma        lsit   etime
 2345   +1.234e-03   +4.321e-03   +5.678e-04   +3.456e+00   +3.456e+00   +1.234e-02   1234   1.23e-02
  SUBit pinfeas      dinfeas      relgap       pobj         dobj         sigma        lsit   etime
   2345 -1.200e-03   -4.100e-03   -5.800e-04   -3.600e+00   -3.600e+00   +1.400e-02      1   1.00e-02

 ┌─────┬───────────┬────────────┬────────────┬────────────┬────────────┬────────────┬──────┬──────────┐
 │SUBit│pinfeas    │ dinfeas    │ relgap     │ pobj       │ dobj       │ sigma      │ lsit │ etime    │
 └─────┴───────────┴────────────┴────────────┴────────────┴────────────┴────────────┴──────┴──────────┘
 └ good termination in subproblem: dualinfes = 4.83e-03, gradLyxi = 9.27e-06, tolsub = 2.28e-03       ┘





    ┌────────────────────────────────────────────────────────────────────────────┐
│iter │  pinfeas  dinfeas  relgap  │    pobj          dobj    │ time │ sigma │
=#
