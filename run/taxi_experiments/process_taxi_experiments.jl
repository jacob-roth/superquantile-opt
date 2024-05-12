
function taxi_path_table(res_alm_grid::Dict)
  res = Dict()
  ks = collect(keys(res_alm_grid))
  ks_prob = sort(ks[[k .!= :walltime for k in ks]])
  res[:times] = [res_alm_grid[k][:walltime] for k in ks_prob]
  res[:iters] = [res_alm_grid[k][:iter] for k in ks_prob]
  res[:pinfeas] = [res_alm_grid[k][:pinfeas] for k in ks_prob]
  res[:dinfeas] = [res_alm_grid[k][:dinfeas] for k in ks_prob]
  res[:pobj] = [res_alm_grid[k][:pobj] for k in ks_prob]
  res[:relgap] = [res_alm_grid[k][:relgap] for k in ks_prob]
  res[:k] = [res_alm_grid[k][:k] for k in ks_prob]
  res[:ctime_vphi] = [sum(res_alm_grid[k][:ctime_vphi]) for k in ks_prob]
  res[:ctime_gphi] = [sum(res_alm_grid[k][:ctime_gphi]) for k in ks_prob]
  res[:ctime_Hphi] = [sum(res_alm_grid[k][:ctime_Hphi]) for k in ks_prob]
  res[:ctime_nd] = [sum(res_alm_grid[k][:ctime_nd]) for k in ks_prob]
  res[:ctime_sort] = [sum(res_alm_grid[k][:ctime_sort]) for k in ks_prob]
  res[:ctime_proj] = [sum(res_alm_grid[k][:ctime_proj]) for k in ks_prob]
  res[:ctime_Ax] = [sum(res_alm_grid[k][:ctime_Ax]) for k in ks_prob]
  res[:ctime_Atlam] = [sum(res_alm_grid[k][:ctime_Atlam]) for k in ks_prob]
  res[:ctime_diagnostics] = [sum(res_alm_grid[k][:ctime_diagnostics]) for k in ks_prob]
  res[:ctime_assign] = [sum(res_alm_grid[k][:ctime_assign]) for k in ks_prob]
  res[:err] = [maximum(r) for r in eachrow(hcat(res[:pinfeas], res[:dinfeas], res[:relgap]))]
  return res
end
scientific_notation(x, d::Integer=2) = begin
  if !isnan(x)
    fmt = Printf.Format("%0.$(d)e")
    s = Printf.format(fmt, x)
    base, expo = split(s, "e")
    out = base * "e" * string(parse(Int64, expo))
  else
    out = L"$-$"
  end
  return out
end

#
# warm
#

outdict = taxi_path_table(res_alm_grid)
nexp = 22
cat = [L"k = \tau\%m", "time [s]", L"$\eta$", "iter", L"$\hat\nabla^2\varphi$", L"$\nabla\varphi$", "sort", "proj"]
ncat = length(cat)
function table_line(i::Integer, res::Dict, pctgrid::Vector)
  line = [
    scientific_notation(pctgrid[i], 2),
    scientific_notation(res[:times][i], 1),
    scientific_notation(res[:err][i], 1),
    @sprintf("%d", res[:iters][i]),
    @sprintf("%0.0f", 100(res[:ctime_Hphi][i] + res[:ctime_nd][i]) / res[:times][i]),
    @sprintf("%0.0f", 100(res[:ctime_gphi][i]) / res[:times][i]),
    @sprintf("%0.0f", 100(res[:ctime_sort][i]) / res[:times][i]),
    @sprintf("%0.0f", 100(res[:ctime_proj][i]) / res[:times][i]),
  ]
end
function total_line(idx, res::Dict, pctgrid::Vector)
  line = [
    L"--",
    scientific_notation(sum(res[:times][idx]), 1),
    L"--",
    L"--",
    @sprintf("%0.0f", 100(sum(res[:ctime_Hphi][idx] + res[:ctime_nd][idx])) / sum(res[:times][idx])),
    @sprintf("%0.0f", 100(sum(res[:ctime_gphi][idx])) / sum(res[:times][idx])),
    @sprintf("%0.0f", 100(sum(res[:ctime_sort][idx])) / sum(res[:times][idx])),
    @sprintf("%0.0f", 100(sum(res[:ctime_proj][idx])) / sum(res[:times][idx])),
  ]
end
using LaTeXTabulars
using LaTeXStrings               # not a dependency, but works nicely
latex_tabular(RESULTPATH * "/warm_table.tex",
  Tabular("l"^(1+ncat)),
  [Rule(:top),
    cat,
    Rule(:mid),
    [table_line(i, outdict, pctgrid) for i in 1:nexp]...,
    Rule(:mid),
    total_line(1:nexp, outdict, pctgrid),
    Rule(:bottom)
  ]
)

#
# cold
#

function table_line(pct::Integer, m::Integer=32534601, n::Integer=19)
  mtds = ["alm", "grb", "grb_oa", "grb_qr"]
  times = []
  for mtd in mtds
    p = RESULTPATH * "cold/obj_1/k_$(pct)pct__L_1/$(m)_$(n)_1__$(mtd)__walltime_i1.csv"
    if ispath(p)
      v = readdlm(p)
      push!(times, v[1])
    else
      push!(times, L"$-$")
    end
  end
  pinfeas = []
  for mtd in mtds
    p = RESULTPATH * "cold/obj_1/k_$(pct)pct__L_1/$(m)_$(n)_1__$(mtd)__pinfeas.csv"
    if ispath(p)
      v = readdlm(p)
      push!(pinfeas, v[1])
    else
      push!(pinfeas, L"$-$")
    end
  end
  dinfeas = []
  for mtd in mtds
    p = RESULTPATH * "cold/obj_1/k_$(pct)pct__L_1/$(m)_$(n)_1__$(mtd)__dinfeas.csv"
    if ispath(p)
      v = readdlm(p)
      push!(dinfeas, v[1])
    else
      push!(dinfeas, L"$-$")
    end
  end
  relgap = []
  for mtd in mtds
    p = RESULTPATH * "cold/obj_1/k_$(pct)pct__L_1/$(m)_$(n)_1__$(mtd)__relgap.csv"
    if ispath(p)
      v = readdlm(p)
      push!(relgap, v[1])
    else
      push!(relgap, L"$-$")
    end
  end
  errs = [max(pinfeas[i], dinfeas[i], relgap[i]) for i in eachindex(mtds)]
  line = [@sprintf("%d", pct/10), "", scientific_notation.(times, 1)..., "", scientific_notation.(errs, 1)...]
  return line
end

bigcat = [L"k = \tau\%m", "", MultiColumn(4, :c, "time [s]"), MultiColumn(4, :c, L"$\eta$")]
smallcat = ["", "", ["p-ALM", "GRB", "GRB-OA", "GRB-QR"]..., "", ["p-ALM", "GRB", "GRB-OA", "GRB-QR"]...]
ncat = length(smallcat)
latex_tabular(RESULTPATH * "/cold_table.tex",
  Tabular("l"^(1+ncat)),
  [Rule(:top),
    bigcat,
    CMidRule(3, 6), CMidRule(8, 11), 
    smallcat,
    Rule(:mid),
    [table_line(10, 32534601, 19), table_line(100, 32534601, 19), table_line(500, 32534601, 19)]...,
    Rule(:mid),
    Rule(:bottom)
  ]
)
