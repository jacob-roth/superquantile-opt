import Pkg; Pkg.activate("."); Pkg.instantiate()
using LinearAlgebra, Random, DelimitedFiles, LaTeXStrings, LaTeXTabulars, Printf, NaNStatistics
using JuMP
include("../../src/SSNALM.jl")
include("helper_alm_experiments.jl")
include("plot_alm_experiments_NEW.jl")
fig, ax = subplots(figsize=(5,5), layout="constrained")
rc("font", size=16)

#
# process setup
#

# global const DATAPATH = "/Users/jakeroth/Downloads/alm_results_20231208/"
# global const DATAPATH = "/Users/jakeroth/Downloads/alm_results_20231213/"
# global const DATAPATH = "/Users/jakeroth/Downloads/alm_results_20240105/"
# global const MAX_N_SHOW = 1_000
# global const MIN_M_SHOW = 100_000
global const DATAPATH = "/Users/jakeroth/Downloads/alm_results_20240131_FULL4/"
global const DATAPATH = "/Users/jakeroth/Downloads/alm_results_5_20240410/"
global const MAX_N_SHOW = 10_000_000
global const MIN_M_SHOW = 0
global const base = 2

#
# formatted printing
#

print_sci_10_digits(d,digits=1) = begin
  """one-exponent-digit and one-base-decimaldigit in base 10"""
  # expo = exponent(d)
  # base = significand(d)
  # sgn = sign(d)
  # @assert(d == sgn*abs(base)*2.0^expo)
  # base10 = 
  # expo10 = 
  # expo = log10(abs(d))
  # base = d * 10^-expo
  # @sprintf("%0.1fe%d", base, expo)
  if d==0
    L"0"
  elseif isnan(d) || isinf(d)
    L"-"
  elseif d <= 100.0 && d >= 0.1
    if d <= 10.0 && d >= 0.1
      if digits==0
        @sprintf("%0.0f", d)
      elseif digits==1
        @sprintf("%0.1f", d)
      elseif digits==2
        @sprintf("%0.1f", d)
      end
    else
      if digits==0
        @sprintf("%0.0f", d)
      elseif digits==1
        @sprintf("%0.1f", d)
      elseif digits==2
        @sprintf("%0.1f", d)
      end
    end
  elseif digits==0
    st = @sprintf("%0.0e", d)
    base, expo = split(st, "e")
    expo = parse(Int64, expo)
    return "$(base)e$(expo)"
  elseif digits==1
    st = @sprintf("%0.1e", d)
    base, expo = split(st, "e")
    expo = parse(Int64, expo)
    return "$(base)e$(expo)"
  elseif digits==2
    st = @sprintf("%0.2e", d)
    base, expo = split(st, "e")
    expo = parse(Int64, expo)
    return "$(base)e$(expo)"
  end
end

# get_mn_grid(outall[oe][on], p*on, md, "alm", "walltime")
# out = outall[oe][on]
# outname = p*on
# method1 = "md"
# method2 = "alm"
# field = "walltime"
# 1000_100
# `find ./ -name "*.orig" | xargs rm -r` or print names without `rm -r`
function get_mn_grid(out::Dict, outname::String, method1::String, method2::String, field::String, min_m_show::Integer=MIN_M_SHOW, max_n_show::Integer=MAX_N_SHOW)
  println("outname=\"$(outname)\"; method1=\"$(method1)\"; method2=\"$(method2)\"; field=\"$(field)\";")
  fs = readdir(outname*"/")
  fs = fs[[f[1] .!= '.' for f in fs]]
  fs_split = split.(fs, "__")
  mnrs = [x[1] for x in fs_split]
  ms_all = [x[1] for x in split.(mnrs, "_")]
  ns_all = [x[2] for x in split.(mnrs, "_")]
  rs_all = [x[3] for x in split.(mnrs, "_")]
  mns = sort(unique([parse.(Int64, [x[1],x[2]]) for x in split.(mnrs, "_")]))
  ms = sort(parse.(Int64, unique(ms_all)))
  ns = sort(parse.(Int64, unique(ns_all)))
  rs = parse.(Int64, unique(rs_all))
  ms = ms[ms.>=min_m_show]
  ns = ns[ns.<=max_n_show]
  num_m = length(ms)
  num_n = length(ns)
  num_r = length(rs)
  println("(num_m,num_n) = $((num_m,num_n))")
  @assert(num_m == num_n)
  data = zeros(num_m, num_n)
  idx = zeros(Bool, num_m, num_n)
  idx_ok1 = zeros(Bool, num_m, num_n)
  idx_ok2 = zeros(Bool, num_m, num_n)
  if method1 == "grb_oa" || method1 == "grb" || method1 == "grbds" || method1 == "grb_best" || method1 == "psgVAN" || method1 == "psgHELI"
    tol = 1e-8
  else
    tol = 1e-3
  end
  myidx = 0
  for i in 1:num_m
    for j in 1:num_n
      if i+j <= num_n+1
        myidx += 1
        mn_key = "$(ms[i])_$(ns[j])"
        idx[i,j] = true

        # grb_opt = all(out["grb"][mn_key]["pinfeas"] .<= tol) && all(out["grb"][mn_key]["dinfeas"] .<= tol) && all(out["grb"][mn_key]["relgap"] .<= tol)
        # grb_opt = all(out["grb"][mn_key]["pinfeas"] .<= 10tol) && all(out["grb"][mn_key]["dinfeas"] .<= 10tol) && all(out["grb"][mn_key]["relgap"] .<= 10tol) #! one case for GRB slightly off: k_1pct__L_10, i=2, j=3
        grb_opt = all(out["grb"][mn_key]["pinfeas"] .<= tol) && all(out["grb"][mn_key]["dinfeas"] .<= tol) && all(out["grb"][mn_key]["relgap"] .<= tol) #! one case for GRB slightly off: k_1pct__L_10, i=2, j=3
        grbds_opt = all(out["grbds"][mn_key]["pinfeas"] .<= tol) && all(out["grbds"][mn_key]["dinfeas"] .<= tol) && all(out["grbds"][mn_key]["relgap"] .<= tol)  
        if method1 == "psgVAN" || method1 == "psgHELI"
          # idx_ok1[i,j] = out[method1][mn_key]["status"][1] == "optimal" # or could be "feasible"
          idx_ok1[i,j] = out[method1][mn_key]["pinfeas"][1] <= 1e-8 && out[method1][mn_key]["pobj"][1] <= out["alm"][mn_key]["pobj"][1] + 1e-4 # or could be "feasible"
        elseif method1 == "grb_best"
          idx_ok1[i,j] = (grb_opt || grbds_opt)
        else
          # idx_ok1[i,j] = out[method1][mn_key]["status"][1] == 1
          idx_ok1[i,j] = all(out[method1][mn_key]["pinfeas"] .<= tol) && all(out[method1][mn_key]["dinfeas"] .<= tol) && all(out[method1][mn_key]["relgap"] .<= tol)
        end

        if method1 == "grb_best"
          if grb_opt && grbds_opt
            data[i,j] = nanmean(min.(out["grb"][mn_key][field], out["grbds"][mn_key][field]) ./ out[method2][mn_key][field])
          elseif grb_opt && !grbds_opt
            data[i,j] = nanmean(out["grb"][mn_key][field] ./ out[method2][mn_key][field])
          elseif !grb_opt && grbds_opt
            data[i,j] = nanmean(out["grbds"][mn_key][field] ./ out[method2][mn_key][field])
          else
            data[i,j] = nanmean(min.(out["grb"][mn_key][field], out["grbds"][mn_key][field]) ./ out[method2][mn_key][field])
          end
        elseif method2 == "grb_best"
          if grb_opt && grbds_opt
            data[i,j] = nanmean(out[method1][mn_key][field] ./ min.(out["grb"][mn_key][field], out["grbds"][mn_key][field]))
          elseif grb_opt && !grbds_opt
            data[i,j] = nanmean(out[method1][mn_key][field] ./ out["grb"][mn_key][field])
          elseif !grb_opt && grbds_opt
            data[i,j] = nanmean(out[method1][mn_key][field] ./ out["grbds"][mn_key][field])
          else
            data[i,j] = nanmean(out[method1][mn_key][field] ./ min.(out["grb"][mn_key][field], out["grbds"][mn_key][field]))
          end
        else
          data[i,j] = nanmean(out[method1][mn_key][field] ./ out[method2][mn_key][field])
        end

        println("(i,j) = $((i,j))")
        println(data[i,j])
        
        if method2 == "psgVAN" || method2 == "psgHELI"
          # idx_ok2[i,j] = out[method2][mn_key]["status"][1] == "optimal" # or could be "feasible"
          idx_ok2[i,j] = out[method2][mn_key]["pinfeas"][1] <= 1e-8 && out[method2][mn_key]["pobj"][1] <= out["alm"][mn_key]["pobj"][1] + 1e-4 # or could be "feasible"
        elseif method2 == "grb_best"
          println("method1 = $method1")
          println("method2 = $method2")
          idx_ok2[i,j] = (grb_opt || grbds_opt)
        else
          # idx_ok2[i,j] = out[method2][mn_key]["status"][1] == 1
          idx_ok2[i,j] = all(out[method2][mn_key]["pinfeas"] .<= tol) && all(out[method2][mn_key]["dinfeas"] .<= tol) && all(out[method2][mn_key]["relgap"] .<= tol)
        end
      end
    end
  end
  return (name=outname, data=data, idx=idx, idx_ok1=idx_ok1, idx_ok2=idx_ok2, ms=ms, ns=ns)
  # idx_ok1 = numerator
  # idx_ok2 = denominator
end

#
# load results
#

mtd = ["alm", "grb_best", "grb", "grbds", "grb_oa", "psgVAN", "psgHELI", "osqp"]
mtd_name = ["ALM", "GRB", "GRB-BR", "GRB-DS", "G-OA", "PSG", "PSG-HELI", "OSQP"]
outall = Dict()
objexpers = readdir(DATAPATH)
objexpers = objexpers[objexpers .!= ".DS_Store"]
objexpers = objexpers[objexpers .!= "._.DS_Store"]
for md in mtd
  if md == "alm"
    continue
  else
    outall[md] = Dict()
    outall[md] = Dict()
    outall[md]["grid_hi_1_ovr"] = -Inf
    outall[md]["grid_hi_1_sol"] = -Inf
    outall[md]["grid_hi_2_ovr"] = -Inf
    outall[md]["grid_hi_2_sol"] = -Inf
    outall[md]["grid_lo_1_ovr"] = +Inf
    outall[md]["grid_lo_1_sol"] = +Inf
    outall[md]["grid_lo_2_ovr"] = +Inf
    outall[md]["grid_lo_2_sol"] = +Inf
  end
end
for oe in objexpers
  objtype = parse(Int64, match(r"obj_(\d+)", oe).captures[1])
  p = DATAPATH * oe * "/"
  outnames = readdir(p)
  outnames = outnames[outnames .!= ".DS_Store"]
  outnames = outnames[outnames .!= "._.DS_Store"]
  outall[oe] = Dict()
  for on in outnames
    L = parse(Int64, match(r"L_(\d+)", on).captures[1])
    k =  parse(Int64, match(r"k_(\d+)", on).captures[1])
    outall[oe][on] = load_alm_results(p, on)
    for (md,mdn) in zip(mtd, mtd_name)
      if md == "alm"
        continue
      else
        println("objtype=$objtype; on=\"$on\"; md=\"$md\"; mdn=\"$mdn\"")
        min_m_show = Int(ceil(MIN_M_SHOW/L))
        max_n_show = Int(MAX_N_SHOW)

        outall[oe][on]["walltime_$(md)_alm"] = get_mn_grid(outall[oe][on], p*on, md, "alm", "walltime", min_m_show, max_n_show)
        outall[oe][on]["walltime_alm_$(md)"] = get_mn_grid(outall[oe][on], p*on, "alm", md, "walltime", min_m_show, max_n_show)
        ratios = copy(outall[oe][on]["walltime_$(md)_alm"][2])
        lose_idx = findall((ratios.<1) .&& (ratios.>0))
        for li in lose_idx
          ratios[li] = -1.0 / ratios[li]
        end
        idx1 = outall[oe][on]["walltime_$(md)_alm"].idx_ok1 # numerator
        idx2 = outall[oe][on]["walltime_$(md)_alm"].idx_ok2 # denominator
        outall[oe][on]["walltime_$(md)_alm_grid"] = ratios
        outall[oe][on]["walltime_alm_idx"] = idx2
        outall[oe][on]["walltime_$(md)_idx"] = idx1
        if objtype == 1
          outall[md]["grid_hi_1_ovr"] = max(outall[md]["grid_hi_1_ovr"], maximum(outall[oe][on]["walltime_$(md)_alm"][2], init=-Inf))
          outall[md]["grid_hi_1_sol"] = max(outall[md]["grid_hi_1_sol"], maximum(outall[oe][on]["walltime_$(md)_alm"].data[idx1 .&& idx2], init=-Inf))
          outall[md]["grid_lo_1_ovr"] = min(outall[md]["grid_lo_1_ovr"], -maximum(outall[oe][on]["walltime_alm_$(md)"][2], init=+Inf))
          outall[md]["grid_lo_1_sol"] = min(outall[md]["grid_lo_1_sol"], -maximum(outall[oe][on]["walltime_alm_$(md)"].data[idx1 .&& idx2], init=+Inf))
        elseif objtype == 2
          outall[md]["grid_hi_2_ovr"] = max(outall[md]["grid_hi_2_ovr"], maximum(outall[oe][on]["walltime_$(md)_alm"][2], init=-Inf))
          outall[md]["grid_hi_2_sol"] = max(outall[md]["grid_hi_2_sol"], maximum(outall[oe][on]["walltime_$(md)_alm"].data[idx1 .&& idx2], init=-Inf))
          outall[md]["grid_lo_2_ovr"] = min(outall[md]["grid_lo_2_ovr"], -maximum(outall[oe][on]["walltime_alm_$(md)"][2], init=+Inf))
          outall[md]["grid_lo_2_sol"] = min(outall[md]["grid_lo_2_sol"], -maximum(outall[oe][on]["walltime_alm_$(md)"].data[idx1 .&& idx2], init=+Inf))
        end
      end
    end
  end
end


# frac optimal PSG: VAN/CAR
opt1 = [outall["obj_1"][k1]["psgVAN"][k2]["status_internal"][1]=="optimal" for k1 in keys(outall["obj_1"]) for k2 in keys(outall["obj_1"][k1]["psgVAN"])];
opt2 = [outall["obj_2"][k1]["psgVAN"][k2]["status_internal"][1]=="optimal" for k1 in keys(outall["obj_2"]) for k2 in keys(outall["obj_2"][k1]["psgVAN"])];
sum(opt1)/length(opt1)
sum(opt2)/length(opt2)

# frac "optimal" PSG: HELI/HELI, not necc optimal numerically...
opt11 = [outall["obj_1"][k1]["psgHELI"][k2]["status_internal"][1]=="optimal" for k1 in keys(outall["obj_1"]) for k2 in keys(outall["obj_1"][k1]["psgHELI"])];
opt22 = [outall["obj_2"][k1]["psgHELI"][k2]["status_internal"][1]=="optimal" for k1 in keys(outall["obj_2"]) for k2 in keys(outall["obj_2"][k1]["psgHELI"])];
sum(opt11)/length(opt11)
sum(opt22)/length(opt22)

#
# plot grid
#

for oe in objexpers
  objtype = parse(Int64, match(r"obj_(\d+)", oe).captures[1])
  p = DATAPATH * oe * "/"
  outnames = readdir(p)
  outnames = outnames[outnames .!= ".DS_Store"]
  outnames = outnames[outnames .!= "._.DS_Store"]
  for on in outnames
    # plot grids
    L = parse(Int64, match(r"L_(\d+)", on).captures[1])
    k =  parse(Int64, match(r"k_(\d+)", on).captures[1])
    for (md,mdn) in zip(mtd, mtd_name)
      if md == "alm"
        continue
      else
        plottitle = "$(mdn)/p-ALM: " * L"L=%$(L)" * ", " *  L"k=%$(k)\%m"
        plotname = "mnheatmap__$(on)__$(objtype)__$(mdn)"
        # vhi = ceil(max(maximum(outall[oe][on]["walltime_$(md)_alm_grid"][end-3:end,1:3]), 1.0))
        # vlo = floor(min(minimum(outall[oe][on]["walltime_$(md)_alm_grid"]), -1.0))
        N = size(outall[oe][on]["walltime_$(md)_alm_grid"],1)
        IDX_i = [N, N-1, N-2, N-3]
        IDX_j = [1, 2, 3, 4]
        vhi = ceil(max(maximum([outall[oe][on]["walltime_$(md)_alm_grid"][i,j] for (i,j) in zip(IDX_i, IDX_j)]), 1.0))
        vlo = floor(min(minimum(outall[oe][on]["walltime_$(md)_alm_grid"]), -1.0))
        plot_grid2(
          # ax,
          outall[oe][on]["walltime_$(md)_alm_grid"],
          outall[oe][on]["walltime_alm_$(md)"].idx,
          outall[oe][on]["walltime_alm_idx"],
          outall[oe][on]["walltime_$(md)_idx"],
          plottitle, plotname, "coolwarm", vlo, vhi,
          outall[oe][on]["walltime_$(md)_alm"].ms,
          outall[oe][on]["walltime_$(md)_alm"].ns,
          base
        )
      end
    end
  end
end

#
# get tables
#

# print_digit(x) = @sprintf("%d",x)
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
function table_compare_line(objtype::Int64, m::Int64, n::Int64, L::Int64, k::Int64, mtds::AbstractVector, alminfo::Bool=false, nrep::Integer=1, repi::Integer=1, base::Integer=2)
  mn = "$(m)_$(n)"
  outname = "k_$(k)pct__L_$(L)/"
  p = DATAPATH * "obj_$(objtype)" * "/" * outname
  
  alm_t_repi = [sum(readdlm(p*"/"*mn*"_$(repi)__alm__walltime.csv")) for repi in 1:nrep]
  alm_t = mean(alm_t_repi)
  alm_t_sort_pct = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__time_sort.csv")) for repi in 1:nrep] ./ alm_t_repi)
  # alm_t_dirn_pct = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__time_dirn.csv")) + sum(readdlm(p*"/"*mn*"_$(repi)__alm__time_Hphi.csv")) for repi in 1:nrep] ./ alm_t_repi)
  alm_t_dirn_pct = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__time_dirn.csv")) for repi in 1:nrep] ./ alm_t_repi)
  alm_t_grad_pct = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__time_gphi.csv")) for repi in 1:nrep] ./ alm_t_repi)
  alm_beta = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__avg_nd_size.csv")) for repi in 1:nrep])
  pobj = [isfile(p*"/"*mn*"_$(repi)__$(mtdi)__pobj.csv") ? nanmean([sum(readdlm(p*"/"*mn*"_$(repi)__$(mtdi)__pobj.csv")) for repi in 1:nrep]) : NaN for mtdi in mtds]
  # pobj_pct = (pobj ./ pobj[mtds.=="alm"]) .- 1
  # pobj_pct = ((pobj ./ abs(pobj[mtds.=="alm"])) .- 1)[mtds.!="alm"]
  pobj_pct = ((pobj .- pobj[mtds.=="alm"]) ./ abs.(pobj[mtds.=="alm"]))[mtds.!="alm"]
  pinfeas = [isfile(p*"/"*mn*"_$(repi)__$(mtdi)__pinfeas.csv") ? nanmean([sum(readdlm(p*"/"*mn*"_$(repi)__$(mtdi)__pinfeas.csv")) for repi in 1:nrep]) : NaN for mtdi in mtds]
  dinfeas = [isfile(p*"/"*mn*"_$(repi)__$(mtdi)__dinfeas.csv") ? nanmean([sum(readdlm(p*"/"*mn*"_$(repi)__$(mtdi)__dinfeas.csv")) for repi in 1:nrep]) : NaN for mtdi in mtds[.!occursin.("psg", mtds)]]
  relgap = [isfile(p*"/"*mn*"_$(repi)__$(mtdi)__relgap.csv") ? nanmean([sum(readdlm(p*"/"*mn*"_$(repi)__$(mtdi)__relgap.csv")) for repi in 1:nrep]) : NaN for mtdi in mtds[.!occursin.("psg", mtds)]]
  walltime = [isfile(p*"/"*mn*"_$(repi)__$(mtdi)__walltime.csv") ? nanmean([sum(readdlm(p*"/"*mn*"_$(repi)__$(mtdi)__walltime.csv")) for repi in 1:nrep]) : NaN for mtdi in mtds]
  linestats = [
    # problem
    base == 2 ? L"(2^{%$(Int64(log2(m)))}, 2^{%$(Int64(log2(n)))})" : L"(10^{%$(Int64(log10(m)))}, 10^{%$(Int64(log10(n)))})",
    "",
    # pobj
    [scientific_notation(x, 1) for x in pobj_pct]...,
    "",
    # pinfeas
    [scientific_notation(x, 1) for x in pinfeas]...,
    "",
    # dinfeas: all except PSG
    [scientific_notation(x, 1) for x in dinfeas]...,
    "",
    # relgap: all except PSG
    [scientific_notation(x, 1) for x in relgap]...,
    "",
    # time
    [scientific_notation(x, 1) for x in walltime]...,
    "",
    # ALM info
    "$(print_sci_10_digits(100alm_t_sort_pct,0))" * L"$\,\lvert\,$" * "$(print_sci_10_digits(100alm_t_dirn_pct,0))" * L"$\,\lvert\,$" * "$(print_sci_10_digits(100alm_t_grad_pct,0))" * L"$\,\lvert\,$" * "$(print_sci_10_digits(alm_beta))"
  ][1:end-2*(1-alminfo)]
  breaks = findall("" .== linestats)
  return linestats, breaks
end

function table_compare_lines(objtype::Int64, L::Int64, k::Int64, mtds::AbstractVector, idx::Union{Nothing,AbstractVector}=nothing, alminfo::Bool=false, nrep::Integer=1, repi::Integer=1)
  outname = "k_$(k)pct__L_$(L)/"
  p = DATAPATH * "obj_$(objtype)" * "/" * outname
  fs = readdir(p)
  fs_split = split.(fs, "__")
  mnrs = [x[1] for x in fs_split]
  ms_all = [x[1] for x in split.(mnrs, "_")] # m
  ns_all = [x[2] for x in split.(mnrs, "_")] # n
  rs_all = [x[3] for x in split.(mnrs, "_")] # rep
  mns = sort(unique([parse.(Int64, [x[1],x[2]]) for x in split.(mnrs, "_")]), rev=true)
  sigm = sortperm([x[1] for x in mns], rev=true)
  sign = sortperm([x[2] for x in mns[sigm]])
  mns = mns[sigm][sign]
  max_n_show = MAX_N_SHOW
  min_m_show = Int(MIN_M_SHOW/L)
  mns = mns[[(x[1] >= min_m_show && x[2] <= max_n_show) for x in mns]]
  if !isnothing(idx)
    # idx = findall(prod.(mns) .== maximum(prod.(mns)))
    mns = mns[idx]
  end
  return [table_compare_line(objtype, mn[1], mn[2], L, k, mtds, alminfo, nrep, repi)[1] for mn in mns]
end

mtd = ["alm", "grb", "grb_oa", "psgHELI", "osqp"]
mtd_name = ["ALM", "GRB", "G-OA", "PSG", "OSQP"]

objtype = 2
m = 1048576
n = 128
L = 1
k = 1
_, breaks = table_compare_line(objtype, m, n, L, k, mtd)
breaks = [0; breaks]
# breaks = [0; 2; 7; 13; 18; 23; 29; 30]
breaks = [0; 2; 7; 13; 18; 23; 28]
mtddict = Dict()
for (k,v) in zip(mtd, mtd_name)
  mtddict[k] = v
end

alminfo = false
idx = [1, 4, 19] # triangle with top m
# idx = [1, 8, 14, 19, 23, 26, 28] # largest
for oe in objexpers
  objtype = parse(Int64, match(r"obj_(\d+)", oe).captures[1])
  p = DATAPATH * oe * "/"
  outnames = readdir(p)
  outnames = outnames[outnames .!= ".DS_Store"]
  outnames = outnames[outnames .!= "._.DS_Store"]
  println("p=$p, outnames=$outnames")
  Lks = [(parse(Int64, match(r"L_(\d+)", on).captures[1]), parse(Int64, match(r"k_(\d+)", on).captures[1])) for on in outnames]
  sig_L = sortperm([x[1] for x in Lks])
  sig_k = sortperm(Lks[sig_L])
  Lks = Lks[sig_L][sig_k]
  latex_tabular(
    PROJPATH*"figures/alm_tables/table_compare__$(oe).tex",
    # Tabular("l"^(breaks[end-(alminfo ? 0 : 1)]-1)*(alminfo ? "c" : "")),
    Tabular("l"^28),
    [
      Rule(:top),
      [
        MultiColumn(1, :c, "Problem"), "",
        MultiColumn(length(mtd)-1, :c, "Relative Pobj"), "",
        MultiColumn(length(mtd), :c, "Pinfeas"), "",
        MultiColumn(length(mtd)-1, :c, "Dinfeas"), "",
        MultiColumn(length(mtd)-1, :c, "Relgap"), "",
        MultiColumn(length(mtd), :c, "Walltime [s]")
        # , "",
        # MultiColumn(1, :c, "ALM Info"), "",
        # alminfo ? MultiColumn(1, :c, "ALM Info") : "",
      ][1:end],
      vcat(
        [
          CMidRule(nothing, nothing, 1, 1),
          CMidRule(nothing, nothing, 3, 6),
          CMidRule(nothing, nothing, 8, 12),
          CMidRule(nothing, nothing, 14, 17),
          CMidRule(nothing, nothing, 19, 22),
          CMidRule(nothing, nothing, 24, 28),
        ][1:end]
      )...,
      [
        L"($m,n$)", "", # problem
        [mtddict[mtdi] for mtdi in mtd[mtd.!="alm"]]..., "", # pobj
        [mtddict[mtdi] for mtdi in mtd]..., "", # pinfeas
        [mtddict[mtdi] for mtdi in mtd[.!occursin.("psg", mtd)]]..., "", # dinfeas
        [mtddict[mtdi] for mtdi in mtd[.!occursin.("psg", mtd)]]..., "", # relgap
        [mtddict[mtdi] for mtdi in mtd]..., "", # walltime
        "w"*L"$\,\lvert\,$"*"x"*L"$\,\lvert\,$"*"y"*L"$\,\lvert\,$"*"z"
        # alminfo ? "w"*L"$\,\lvert\,$"*"x"*L"$\,\lvert\,$"*"y"*L"$\,\lvert\,$"*"z" : "",  # alm info
      ][1:end-2*(1-alminfo)],
      Rule(:mid),
      vcat([[Rule(:mid), [MultiColumn(28, :c, L"L={%$(Lk[1])},\;k={%$(Lk[2])}\%m")], Rule(:mid), table_compare_lines(objtype, Lk[1], Lk[2], mtd, idx)...] for Lk in Lks]...)...,
      Rule(:bottom),
    ]
  )
end
latex_tabular(
  stdout,
  Tabular("l"^31),
  [CMidRule("0.1em", nothing, 1, 31),""]
)

#
# ALM 2nd-order-info tables combine 1/2
#

function table_detail_line(m::Int64, n::Int64, L::Int64, k::Int64, mtds::AbstractVector, alminfo::Bool=false, nrep::Integer=1, repi::Integer=1, base::Integer=2)
  mn = "$(m)_$(n)"
  outname = "k_$(k)pct__L_$(L)/" 
  p1 = DATAPATH * "obj_1" * "/" * outname
  p2 = DATAPATH * "obj_2" * "/" * outname

  alm_iter__1 = mean([sum(readdlm(p1*"/"*mn*"_$(repi)__alm__iter.csv")) for repi in 1:nrep])
  alm_ssniter__1 = mean([sum(sum(readdlm(p1*"/"*mn*"_$(repi)__alm__nit.csv"))) for repi in 1:nrep])
  alm_t_repi__1 = [sum(readdlm(p1*"/"*mn*"_$(repi)__alm__walltime.csv")) for repi in 1:nrep]
  alm_t__1 = mean(alm_t_repi__1)
  alm_t_sort_pct__1 = mean([sum(readdlm(p1*"/"*mn*"_$(repi)__alm__time_sort.csv")) for repi in 1:nrep] ./ alm_t_repi__1)
  alm_t_dirn_pct__1 = mean([sum(readdlm(p1*"/"*mn*"_$(repi)__alm__time_dirn.csv")) for repi in 1:nrep] ./ alm_t_repi__1)
  # alm_t_dirn_pct__1 = mean([sum(readdlm(p1*"/"*mn*"_$(repi)__alm__time_dirn.csv")) + sum(readdlm(p1*"/"*mn*"_$(repi)__alm__time_Hphi.csv")) for repi in 1:nrep] ./ alm_t_repi__1) #! dirn has both already in output; confusing...
  alm_t_grad_pct__1 = mean([sum(readdlm(p1*"/"*mn*"_$(repi)__alm__time_gphi.csv")) for repi in 1:nrep] ./ alm_t_repi__1)
  alm_beta__1 = sum(mean([nanmean(readdlm(p1*"/"*mn*"_$(repi)__alm__absbeta.csv")) for repi in 1:nrep]) * sum([readdlm(p1*"/"*mn*"_$(repi)__alm__cvar.csv") .>= readdlm(p1*"/"*mn*"_$(repi)__alm__cvar_rhs.csv").-1e-5 for repi in 1:nrep])) # first part is average per active L; multiply by active L
  pinfeas__1 = mean([sum(readdlm(p1*"/"*mn*"_$(repi)__alm__pinfeas.csv")) for repi in 1:nrep])
  dinfeas__1 = mean([sum(readdlm(p1*"/"*mn*"_$(repi)__alm__dinfeas.csv")) for repi in 1:nrep])
  relgap__1 = mean([sum(readdlm(p1*"/"*mn*"_$(repi)__alm__relgap.csv")) for repi in 1:nrep])
  kktres__1 = max(pinfeas__1, dinfeas__1, relgap__1)
  walltime__1 = mean([sum(readdlm(p1*"/"*mn*"_$(repi)__alm__walltime.csv")) for repi in 1:nrep])
  walltime_i1__1 = mean([sum(readdlm(p1*"/"*mn*"_$(repi)__alm__walltime_i1.csv")) for repi in 1:nrep])
  walltime_i2__1 = mean([sum(readdlm(p1*"/"*mn*"_$(repi)__alm__walltime_i2.csv")) for repi in 1:nrep])
  avg_nd_size__1 = mean([mean(readdlm(p1*"/"*mn*"_$(repi)__alm__avg_nd_size.csv")) for repi in 1:nrep])
  alm_final_beta_1 = sum(mean([(readdlm(p1*"/"*mn*"_$(repi)__alm__absbeta.csv"))[end] for repi in 1:nrep]) * sum([readdlm(p1*"/"*mn*"_$(repi)__alm__cvar.csv") .>= readdlm(p1*"/"*mn*"_$(repi)__alm__cvar_rhs.csv").-1e-5 for repi in 1:nrep]))
  
  alm_iter__2 = mean([sum(readdlm(p2*"/"*mn*"_$(repi)__alm__iter.csv")) for repi in 1:nrep])
  alm_ssniter__2 = mean([sum(sum(readdlm(p2*"/"*mn*"_$(repi)__alm__nit.csv"))) for repi in 1:nrep])
  alm_t_repi__2 = [sum(readdlm(p2*"/"*mn*"_$(repi)__alm__walltime.csv")) for repi in 1:nrep]
  alm_t__2 = mean(alm_t_repi__2)
  alm_t_sort_pct__2 = mean([sum(readdlm(p2*"/"*mn*"_$(repi)__alm__time_sort.csv")) for repi in 1:nrep] ./ alm_t_repi__2)
  # alm_t_dirn_pct__2 = mean([sum(readdlm(p2*"/"*mn*"_$(repi)__alm__time_dirn.csv")) + sum(readdlm(p2*"/"*mn*"_$(repi)__alm__time_Hphi.csv")) for repi in 1:nrep] ./ alm_t_repi__2) #! dirn has both already in output; confusing...
  alm_t_dirn_pct__2 = mean([sum(readdlm(p2*"/"*mn*"_$(repi)__alm__time_dirn.csv")) for repi in 1:nrep] ./ alm_t_repi__2)
  alm_t_grad_pct__2 = mean([sum(readdlm(p2*"/"*mn*"_$(repi)__alm__time_gphi.csv")) for repi in 1:nrep] ./ alm_t_repi__2)
  alm_beta__2 = sum(mean([nanmean(readdlm(p2*"/"*mn*"_$(repi)__alm__absbeta.csv")) for repi in 1:nrep]) * sum([readdlm(p2*"/"*mn*"_$(repi)__alm__cvar.csv") .>= readdlm(p2*"/"*mn*"_$(repi)__alm__cvar_rhs.csv").-1e-5 for repi in 1:nrep])) # first part is average per active L; multiply by active L
  pinfeas__2 = mean([sum(readdlm(p2*"/"*mn*"_$(repi)__alm__pinfeas.csv")) for repi in 1:nrep])
  dinfeas__2 = mean([sum(readdlm(p2*"/"*mn*"_$(repi)__alm__dinfeas.csv")) for repi in 1:nrep])
  relgap__2 = mean([sum(readdlm(p2*"/"*mn*"_$(repi)__alm__relgap.csv")) for repi in 1:nrep])
  kktres__2 = max(pinfeas__2, dinfeas__2, relgap__2)
  walltime__2 = mean([sum(readdlm(p2*"/"*mn*"_$(repi)__alm__walltime.csv")) for repi in 1:nrep])
  walltime_i1__2 = mean([sum(readdlm(p2*"/"*mn*"_$(repi)__alm__walltime_i1.csv")) for repi in 1:nrep])
  walltime_i2__2 = mean([sum(readdlm(p2*"/"*mn*"_$(repi)__alm__walltime_i2.csv")) for repi in 1:nrep])
  avg_nd_size__2 = mean([mean(readdlm(p2*"/"*mn*"_$(repi)__alm__avg_nd_size.csv")) for repi in 1:nrep])
  alm_final_beta_2 = sum(mean([(readdlm(p2*"/"*mn*"_$(repi)__alm__absbeta.csv"))[end] for repi in 1:nrep]) * sum([readdlm(p2*"/"*mn*"_$(repi)__alm__cvar.csv") .>= readdlm(p2*"/"*mn*"_$(repi)__alm__cvar_rhs.csv").-1e-5 for repi in 1:nrep]))
  
  linestats = [
    # problem
    base == 2 ? L"(2^{%$(Int64(log2(m)))}, 2^{%$(Int64(log2(n)))})" : L"(10^{%$(Int64(log10(m)))}, 10^{%$(Int64(log10(n)))})",
    "",
    # pinfeas/dinfeas/relgap
    [scientific_notation(x, 1) for x in kktres__1]...,
    [scientific_notation(x, 1) for x in kktres__2]...,
    "",
    # timei1
    [scientific_notation(x, 1) for x in walltime_i1__1]...,
    [scientific_notation(x, 1) for x in walltime_i1__2]...,
    "",
    # timei2
    [scientific_notation(x, 1) for x in walltime_i2__1]...,
    [scientific_notation(x, 1) for x in walltime_i2__2]...,
    "",
    # time
    [scientific_notation(x, 1) for x in walltime__1]...,
    [scientific_notation(x, 1) for x in walltime__2]...,
    "",
    # ALM info
    "$(print_sci_10_digits(100alm_t_sort_pct__1,0))",
    "$(print_sci_10_digits(100alm_t_sort_pct__2,0))",
    "",
    "$(print_sci_10_digits(100alm_t_dirn_pct__1,0))",
    "$(print_sci_10_digits(100alm_t_dirn_pct__2,0))",
    "",
    "$(print_sci_10_digits(100alm_t_grad_pct__1,0))",
    "$(print_sci_10_digits(100alm_t_grad_pct__2,0))",
    # "",
    # (avg_nd_size__1<10000 ? "$(Int(round(avg_nd_size__1,digits=-1)))" : "$(print_sci_10_digits(avg_nd_size__1))"),
    # (avg_nd_size__2<10000 ? "$(Int(round(avg_nd_size__2,digits=-1)))" : "$(print_sci_10_digits(avg_nd_size__2))"),
    "",
    (alm_ssniter__1<95 ? "$(Int(round(alm_iter__1,digits=0))) "*L"\mid"*" $(Int(round(alm_ssniter__1,digits=-1)))" * L"\phantom{0}" : "$(Int(round(alm_iter__1,digits=0))) "*L"\mid"*" $(Int(round(alm_ssniter__1,digits=-1)))"),
    (alm_ssniter__2<95 ? "$(Int(round(alm_iter__2,digits=0))) "*L"\mid"*" $(Int(round(alm_ssniter__2,digits=-1)))" * L"\phantom{0}" : "$(Int(round(alm_iter__2,digits=0))) "*L"\mid"*" $(Int(round(alm_ssniter__2,digits=-1)))"),
    "",
    (alm_final_beta_1 < 10000 ? "$(Int(round(alm_beta__1,digits=-1)))" : "$(print_sci_10_digits(alm_beta__1))"),
    (alm_final_beta_2 < 10000 ? "$(Int(round(alm_beta__2,digits=-1)))" : "$(print_sci_10_digits(alm_beta__2))"),
    "",
    (alm_final_beta_1 < 10000 ? "$(Int(round(alm_final_beta_1,digits=-1)))" : "$(print_sci_10_digits(alm_final_beta_1))"),
    (alm_final_beta_2 < 10000 ? "$(Int(round(alm_final_beta_2,digits=-1)))" : "$(print_sci_10_digits(alm_final_beta_2))"),
  ]
  breaks = findall("" .== linestats)
  return linestats, breaks
end
function table_detail_lines(L::Int64, k::Int64, mtds::AbstractVector, idx::Union{Nothing,AbstractVector}=nothing, alminfo::Bool=false, nrep::Integer=1, repi::Integer=1)
  outname = "k_$(k)pct__L_$(L)/"
  p = DATAPATH * "obj_1" * "/" * outname #! ensure same experiment sizes on 1 and 2 obj
  fs = readdir(p)
  fs_split = split.(fs, "__")
  mnrs = [x[1] for x in fs_split]
  ms_all = [x[1] for x in split.(mnrs, "_")] # m
  ns_all = [x[2] for x in split.(mnrs, "_")] # n
  rs_all = [x[3] for x in split.(mnrs, "_")] # rep
  mns = sort(unique([parse.(Int64, [x[1],x[2]]) for x in split.(mnrs, "_")]), rev=true)
  sigm = sortperm([x[1] for x in mns], rev=true)
  sign = sortperm([x[2] for x in mns[sigm]])
  mns = mns[sigm][sign]
  max_n_show = MAX_N_SHOW
  min_m_show = Int(MIN_M_SHOW/L)
  mns = mns[[(x[1] >= min_m_show && x[2] <= max_n_show) for x in mns]]
  if !isnothing(idx)
    # idx = findall(prod.(mns) .== maximum(prod.(mns)))
    mns = mns[idx]
  end
  return [table_detail_line(mn[1], mn[2], L, k, mtds, alminfo, nrep, repi)[1] for mn in mns]
end

# idx = [1, 4, 19] # triangle with top m
# idx = [1, 8, 14, 19, 23, 26, 28] # largest
# idx = [1, 2, 8] # largest triangle: 2
idx = [1, 2, 3, 8, 9, 14] # largest triangle: 3
idx = [1, 2, 3, 4, 8, 9, 10, 14, 15, 19] # largest triangle: 4
Lks = [
  (1, 1),
  (1, 10),
  (10, 1),
  (10, 10)
]

latex_tabular(
  PROJPATH*"figures/alm_tables/table_detail.tex",
  # Tabular("c"^25),
  Tabular("c"^31),
  [
    Rule(:top),
    [
      MultiColumn(1, :c, "Problem"), "",
      MultiColumn(2, :c, "KKT residual"), "",
      MultiColumn(2, :c, L"$\text{Time}_{\text{a}}: \epsilon=10^{-3}$"), "",
      MultiColumn(2, :c, L"$\text{Time}_{\text{b}}: \epsilon=10^{-6}$"), "",
      MultiColumn(2, :c, L"$\text{Time}_{\text{c}}: \epsilon=10^{-8}$"), "",
      MultiColumn(2, :c, L"\% $\text{Time}_{\text{c}}$: Sort"), "",
      MultiColumn(2, :c, L"\% $\text{Time}_{\text{c}}$: $\widehat{\nabla}^2\varphi\setminus v$"), "",
      MultiColumn(2, :c, L"\% $\text{Time}_{\text{c}}$: $\widehat{\nabla}\varphi$"), "",
      # MultiColumn(2, :c, L"\text{avg SSN dim}"), "",
      MultiColumn(2, :c, L"\text{ALM }\mid\text{ SSN iter}"), "",
      MultiColumn(2, :c, L"\text{avg }|\bar\beta|"), "",
      MultiColumn(2, :c, L"\text{final }|\bar\beta|"),
    ],
    vcat(
      [
        CMidRule(3i, 3i+1)
        for i in 1:10
      ]
    )...,
    vcat(
      L"($m,n$)", "", # problem
      ["lin.", "quad."], "",
      ["lin.", "quad."], "",
      ["lin.", "quad."], "",
      ["lin.", "quad."], "",
      ["lin.", "quad."], "",
      ["lin.", "quad."], "",
      ["lin.", "quad."], "",
      ["lin.", "quad."], "",
      ["lin.", "quad."], "",
      ["lin.", "quad."]
    ),
    Rule(:mid),
    vcat([[Rule(:mid), [MultiColumn(31, :c, L"L={%$(Lk[1])},\;k={%$(Lk[2])}\%m")], Rule(:mid), table_detail_lines(Lk[1], Lk[2], mtd, idx)...] for Lk in Lks]...)...,
    Rule(:bottom),
  ]
)

#=======================
# shell rename commands

# 1 thread
for f in *_alm_mt_*.csv; do mv -v "$f" "${f/alm_mt/alm_mt_1}"; done;
for f in *_grbps_*.csv; do mv -v "$f" "${f/grbps/grbps_1}"; done;
for f in *_grbnps_*.csv; do mv -v "$f" "${f/grbnps/grbnps_1}"; done;

# 4 thread
for f in *_alm_mt_*.csv; do mv -v "$f" "${f/alm_mt/alm_mt_4}"; done;
for f in *_grbps_*.csv; do mv -v "$f" "${f/grbps/grbps_4}"; done;
for f in *_grbnps_*.csv; do mv -v "$f" "${f/grbnps/grbnps_4}"; done;
========================#


#=
#
# ALM 2nd-order-info tables
#


function table_detail_line(objtype::Int64, m::Int64, n::Int64, L::Int64, k::Int64, mtds::AbstractVector, alminfo::Bool=false, nrep::Integer=1, repi::Integer=1, base::Integer=2)
  mn = "$(m)_$(n)"
  outname = "k_$(k)pct__L_$(L)/"
  p = DATAPATH * "obj_$(objtype)" * "/" * outname
  
  alm_t_repi = [sum(readdlm(p*"/"*mn*"_$(repi)__alm__walltime.csv")) for repi in 1:nrep]
  alm_t = mean(alm_t_repi)
  alm_t_sort_pct = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__time_sort.csv")) for repi in 1:nrep] ./ alm_t_repi)
  alm_t_dirn_pct = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__time_dirn.csv")) + sum(readdlm(p*"/"*mn*"_$(repi)__alm__time_Hphi.csv")) for repi in 1:nrep] ./ alm_t_repi)
  alm_t_grad_pct = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__time_gphi.csv")) for repi in 1:nrep] ./ alm_t_repi)
  alm_beta = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__avg_nd_size.csv")) for repi in 1:nrep])
  pinfeas = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__pinfeas.csv")) for repi in 1:nrep])
  dinfeas = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__dinfeas.csv")) for repi in 1:nrep])
  relgap = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__relgap.csv")) for repi in 1:nrep])
  kktres = max(pinfeas, dinfeas, relgap)
  walltime = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__walltime.csv")) for repi in 1:nrep])
  walltime_i1 = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__walltime_i1.csv")) for repi in 1:nrep])
  walltime_i2 = mean([sum(readdlm(p*"/"*mn*"_$(repi)__alm__walltime_i2.csv")) for repi in 1:nrep])
  linestats = [
    # problem
    base == 2 ? L"(2^{%$(Int64(log2(m)))}, 2^{%$(Int64(log2(n)))})" : L"(10^{%$(Int64(log10(m)))}, 10^{%$(Int64(log10(n)))})",
    "",
    # pinfeas/dinfeas/relgap
    [scientific_notation(x, 1) for x in kktres]...,
    "",
    # timei1
    [scientific_notation(x, 1) for x in walltime_i1]...,
    "",
    # timei2
    [scientific_notation(x, 1) for x in walltime_i2]...,
    "",
    # time
    [scientific_notation(x, 1) for x in walltime]...,
    "",
    # ALM info
    "$(print_sci_10_digits(100alm_t_sort_pct,0))",
    "",
    "$(print_sci_10_digits(100alm_t_dirn_pct,0))",
    "",
    "$(print_sci_10_digits(100alm_t_grad_pct,0))",
    "",
    "$(print_sci_10_digits(alm_beta))",
  ]
  breaks = findall("" .== linestats)
  return linestats, breaks
end
function table_detail_lines(objtype::Int64, L::Int64, k::Int64, mtds::AbstractVector, idx::Union{Nothing,AbstractVector}=nothing, alminfo::Bool=false, nrep::Integer=1, repi::Integer=1)
  outname = "k_$(k)pct__L_$(L)/"
  p = DATAPATH * "obj_$(objtype)" * "/" * outname
  fs = readdir(p)
  fs_split = split.(fs, "__")
  mnrs = [x[1] for x in fs_split]
  ms_all = [x[1] for x in split.(mnrs, "_")] # m
  ns_all = [x[2] for x in split.(mnrs, "_")] # n
  rs_all = [x[3] for x in split.(mnrs, "_")] # rep
  mns = sort(unique([parse.(Int64, [x[1],x[2]]) for x in split.(mnrs, "_")]), rev=true)
  sigm = sortperm([x[1] for x in mns], rev=true)
  sign = sortperm([x[2] for x in mns[sigm]])
  mns = mns[sigm][sign]
  max_n_show = MAX_N_SHOW
  min_m_show = Int(MIN_M_SHOW/L)
  mns = mns[[(x[1] >= min_m_show && x[2] <= max_n_show) for x in mns]]
  if !isnothing(idx)
    # idx = findall(prod.(mns) .== maximum(prod.(mns)))
    mns = mns[idx]
  end
  return [table_detail_line(objtype, mn[1], mn[2], L, k, mtds, alminfo, nrep, repi)[1] for mn in mns]
end

# idx = [1, 4, 19] # triangle with top m
# idx = [1, 8, 14, 19, 23, 26, 28] # largest
idx = [1, 2, 8] # largest triangle
# idx = [1, 2, 3, 8, 9, 14] # largest triangle
for oe in objexpers
  objtype = parse(Int64, match(r"obj_(\d+)", oe).captures[1])
  p = DATAPATH * oe * "/"
  outnames = readdir(p)
  outnames = outnames[outnames .!= ".DS_Store"]
  outnames = outnames[outnames .!= "._.DS_Store"]
  println("p=$p, outnames=$outnames")
  Lks = [(parse(Int64, match(r"L_(\d+)", on).captures[1]), parse(Int64, match(r"k_(\d+)", on).captures[1])) for on in outnames]
  sig_L = sortperm([x[1] for x in Lks])
  sig_k = sortperm(Lks[sig_L])
  Lks = Lks[sig_L][sig_k]
  latex_tabular(
    PROJPATH*"figures/alm_tables/table_detail__$(oe).tex",
    Tabular("l"^17),
    [
      Rule(:top),
      [
        L"($m,n$)", "",
        "KKT residual", "",
        L"T: $\epsilon=10^{-3}$", "",
        L"T: $\epsilon=10^{-6}$", "",
        L"T: $\epsilon=10^{-8}$", "",
        "\\%: Sort", "",
        L"\%: $\hat{\nabla}^2\varphi\setminus v$", "",
        L"\%: $\hat{\nabla}\varphi$", "",
        L"|\beta|",
      ],
      vcat(
        [
          CMidRule(2i-1, 2i-1)
          for i in 1:9
        ]
      )...,
      Rule(:mid),
      vcat([[Rule(:mid), [MultiColumn(17, :c, L"L={%$(Lk[1])},\;k={%$(Lk[2])}\%m")], Rule(:mid), table_detail_lines(objtype, Lk[1], Lk[2], mtd, idx)...] for Lk in Lks]...)...,
      Rule(:bottom),
    ]
  )
end
=#