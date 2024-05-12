using PyPlot, PyCall
@pyimport mpl_toolkits.axes_grid1 as ag1
@pyimport matplotlib.patches as mpatches
# @pyimport matplotlib.ticker as ticker
# np = pyimport("numpy")
rc("text", usetex=true)
rc("font", family="serif")
rc("mathtext", fontset="cm")
rc("legend", fancybox=false)
cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

using Printf

function format_and_round(number::Float64)
  # Round the number based on its absolute value
  rounded_number = if abs(number) <= 5
    round(number * 2) / 2
  else
    round(number)
  end

  # Determine formatting based on whether the number is an integer
  if rounded_number == round(rounded_number)
    # Format as an integer
    formatted_string = @sprintf("%d", rounded_number)
  else
    # Format as a decimal, will not exceed 5 characters
    formatted_string = @sprintf("%.1f", rounded_number)
  end

  # Pad the string to ensure it has a length of 5 characters
  padded_string = lpad(formatted_string, 5)

  return padded_string
end

@pydef mutable struct CustomFormatter <: matplotlib.ticker.Formatter
  function __call__(self, x, pos=nothing)
    return format_and_round(x)
  end
end

function plot_grid2(
  # ax,
  ratios::Matrix{Tf}, idx::Matrix{Bool}, idx_alm::Matrix{Bool}, idx_mtd::Matrix{Bool},
  plottitle::String, plotname::String,
  colorscheme::String, vlo::Float64, vhi::Float64,
  ms::Vector{Int64}, ns::Vector{Int64},
  base::Int64=10,
  fs::Integer=20,
) where Tf<:Real
  # get axes
  fig, ax = subplots(figsize=(5,5), layout="constrained")
  rc("font", size=16)
  # ax = subplots(1,1,1)
  ax.cla()
  cmap = get_cmap(colorscheme)
  divnorm = matplotlib.colors.TwoSlopeNorm(vmin=vlo, vcenter=0.0, vmax=vhi)

  im = ax.imshow(clamp.(ratios, vlo, vhi), alpha=1, cmap=cmap, norm=divnorm, origin="lower")
  divider = ag1.axes_divider.make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)

  # colorbar
  cbar = colorbar(im, cax=cax, orientation="vertical", extend="max", format=CustomFormatter())
  cbar[:set_ticks]([floor(vlo), floor(vlo*2)/4, 0, ceil(vhi*2)/4, ceil(vhi)])

  # overwrite inactives
  for i in 1:size(idx, 1)
    for j in 1:size(idx, 2)
      if idx[i, j] == 0
        rect = mpatches.Rectangle((j-1.5, i-1.5), 1, 1, color="white")
        ax.add_patch(rect)
      end
    end
  end
  
  # findall ALM non-optimal
  for ij in findall(.!idx_alm)
    # https://stackoverflow.com/a/14049483/6272122
    i = ij[2]
    j = ij[1]
    println("hatching $ij")
    if idx[ij[1],ij[2]]
      ax.add_patch(matplotlib.patches.Rectangle((i-1.5, j-1.5), 1, 1, hatch="--", linewidth=0, alpha=0.9, fill=false, snap=false, edgecolor=nothing))
    end
  end

  # findall MTD non-optimal
  for ij in findall(.!idx_mtd)
    # https://stackoverflow.com/a/14049483/6272122
    i = ij[2]
    j = ij[1]
    println("hatching $ij")
    if idx[ij[1],ij[2]]
      ax.add_patch(matplotlib.patches.Rectangle((i-1.5, j-1.5), 1, 1, hatch="||", linewidth=0, alpha=0.9, fill=false, snap=false, edgecolor=nothing))
    end
  end

  tlo_m = 0
  thi_m = length(ms)-1
  tlen_m = length(ms)
  tlo_n = 0
  thi_n = length(ns)-1
  tlen_n = length(ns)
  if base == 10
    log_ns = log10.(ns)
    log_ms = log10.(ms)
  elseif base == 2
    log_ns = log2.(ns)
    log_ms = log2.(ms)
  end
  fmt_log_ns=[L"%$(base)^{%$(j)}" for j in Int.(log_ns)]
  fmt_log_ms=[L"%$(base)^{%$(j)}" for j in Int.(log_ms)]
  ax.set_xticks(collect(range(tlo_n, stop=thi_n, length=tlen_n)), fmt_log_ns)
  ax.set_yticks(collect(range(tlo_m, stop=thi_m, length=tlen_m)), fmt_log_ms)
  ax.set_xlabel(L"n", fontsize=fs)
  ax.set_ylabel(L"m_l", fontsize=fs)
  ax.set_frame_on(false)
  ax.set_title(plottitle)
  # https://stackoverflow.com/q/76619274/6272122
  plt.tight_layout()
  plt.subplots_adjust(right = 0.9)
  println(PROJPATH * "figures/alm_heatmap/$(plotname).pdf")
  # fig.savefig(PROJPATH * "figures/alm_heatmap/$(plotname).pdf", bbox_inches = "tight", pad_inches=0.02)
  fig.savefig(PROJPATH * "figures/alm_heatmap/$(plotname).pdf", pad_inches=0.02) # https://stackoverflow.com/q/76619274/6272122
  cbar.remove()
  close(fig)
end

function plot_diagnostics(out::Dict, outname::String, method::String, plotname::String)
  fs = readdir(DATAPATH * outname*"/")
  fs_split = split.(fs, "__")
  pct = parse(Int64, match(r"(\d+)", split(split(outname, "__")[1], "_")[2]).match)/100
  mnrs = [x[1] for x in fs_split]
  ms_all = [x[1] for x in split.(mnrs, "_")]
  ns_all = [x[2] for x in split.(mnrs, "_")]
  rs_all = [x[3] for x in split.(mnrs, "_")]
  mns = sort(unique([parse.(Int64, [x[1],x[2]]) for x in split.(mnrs, "_")]))
  ms = sort(parse.(Int64, unique(ms_all)))
  ns = sort(parse.(Int64, unique(ns_all)))
  rs = parse.(Int64, unique(rs_all))
  num_m = length(ms)
  num_n = length(ns)
  num_r = length(rs)
  @assert(num_m == num_n)
  for n in ns[1:end-1]
    low_n = n
    low_n_mask = findall(x[2] == low_n for x in mns)
    fields = [
      "time_vphi", "time_gphi", "time_dirn",
      "time_sort", "time_proj",#"time_A"#, "time_At",#"walltime"#, "lsit"#, "iter"
    ]
    divs = [
      "lsit", "nit", "nit",
      "lsit", "lsit",# "nit",
    ]
    fieldnames = [
      L"$\varphi$", L"$\nabla\varphi$", L"$\hat{\nabla}^2\varphi\setminus v$",
      "sort", "proj",# L"$Ax$", L"$A^\top \lambda$", "total"# "ALMit"
    ]
    mn_keys = [string(m)*"_"*string(n) for (m,n) in mns[low_n_mask]]
    mn_mask_m = [x[1] for x in mns[low_n_mask]]
    # axx = ax.twinx()
    # for time_type in [:total, :final]
    mindata = 1.0
    for time_type in [:final]
      ax.cla()
      for (f,fn,dv) in zip(fields,fieldnames,divs)
        #! must adjust if have multiple reps
        if time_type == :total
          data = [sum(sum(out[method][mn_key][f])) for mn_key in mn_keys]
          ax.set_title("Cumulative time")
        elseif time_type == :final
          # dvby = ""
          # dvby *= "dividing by: "
          # for mn_key in mn_keys
          #   dvby *= "$dv = $(out[method][mn_key][dv][1][end])\t"
          # end
          # println(dvby)
          # val = ""
          # val *= "val: "
          # for mn_key in mn_keys
          #   val *= "$(out[method][mn_key][f][1][end])\t"
          # end
          # println(val)
          # data = [out[method][mn_key][f][1][end] ./ out[method][mn_key][dv][1][end] for mn_key in mn_keys]
          data = [mean(out[method][mn_key][f][1] ./ out[method][mn_key][dv][1]) for mn_key in mn_keys]
          # data = [out[method][mn_key][f][1][end] ./ out[method][mn_key]["lsit"][1][end] for mn_key in mn_keys]
          # data = [out[method][mn_key][f][1][end] for mn_key in mn_keys]
          # data ./= data[1]
          
          # remove linear trend
          # data ./= (mn_mask_m./mn_mask_m[1])
          
          # minimum
          # mindata = min(minimum(data), mindata)
          # mindata = mean([mean(data[1]), mindata])
          # println(mindata)
          ax.set_title("Average evaluation time")
        end
        ax.loglog(mn_mask_m, data, label=fn, linestyle="--", marker="o")
      end
      # ax.grid(true)
      ax.grid(false)
      ax.set_xlabel(L"$m$")
      leg = ax.legend(loc="upper left")
      # leg = ax.legend(loc="lower left")
      # leg.remove()
      
      # comparisons
      # ax.cla()
      # ax.loglog(mn_mask_m, (mn_mask_m.^2)./(mn_mask_m[1])^2, linestyle="-", label=L"$m^2$", color="black", zorder=1, linewidth=3)
      # axx.loglog(mn_mask_m, mn_mask_m .* log10.(mn_mask_m.*pct), linestyle=":", label=L"$m\log(k)$", color="black", zorder=1, linewidth=3)
      
      # # base trend
      # ax.loglog(mn_mask_m, mindata .* (mn_mask_m .* log10.(mn_mask_m)./(mn_mask_m[1] * log10(mn_mask_m[1]))), linestyle=":", label=L"$m\log(m)$", color="black", zorder=1, linewidth=3)
      # ax.loglog(mn_mask_m, mindata .* (mn_mask_m .* log10.(mn_mask_m.*pct)./(mn_mask_m[1] .* log10(mn_mask_m[1].*pct))), linestyle="-", label=L"$m\log(k)$", color="black", zorder=1, linewidth=3)
      # ax.loglog(mn_mask_m, mindata .* (mn_mask_m./mn_mask_m[1]), linestyle="-", label=L"$m$", color="black", zorder=1, linewidth=3)

      # # remove linear trend
      # ax.loglog(mn_mask_m, mindata .* (mn_mask_m .* log10.(mn_mask_m)./(mn_mask_m[1] * log10(mn_mask_m[1]))) ./ (mn_mask_m./mn_mask_m[1]), linestyle=":", label=L"$m\log(m)$", color="black", zorder=1, linewidth=3)
      # ax.loglog(mn_mask_m, mindata .* (mn_mask_m .* log10.(mn_mask_m.*pct)./(mn_mask_m[1] .* log10(mn_mask_m[1].*pct))) ./ (mn_mask_m./mn_mask_m[1]), linestyle="-", label=L"$m\log(k)$", color="black", zorder=1, linewidth=3)
      # ax.loglog(mn_mask_m, mindata .* (mn_mask_m./mn_mask_m[1]) ./ (mn_mask_m./mn_mask_m[1]), linestyle="-", label=L"$m$", color="black", zorder=1, linewidth=3)

      # comparisons
      # axx.cla()
      # axx.loglog(mn_mask_m, mn_mask_m.^2, linestyle="-", label=L"$m^2$", color="black", zorder=1, linewidth=3)
      # # axx.loglog(mn_mask_m, mn_mask_m .* log10.(mn_mask_m.*pct), linestyle=":", label=L"$m\log(k)$", color="black", zorder=1, linewidth=3)
      # axx.loglog(mn_mask_m, mn_mask_m .* log10.(mn_mask_m), linestyle=":", label=L"$m\log(m)$", color="black", zorder=1, linewidth=3)
      # axx.loglog(mn_mask_m, mn_mask_m, linestyle="-", label=L"$m$", color="black", zorder=1, linewidth=3)
      
      # ax.grid(true)
      # Axtime_blas = zeros(length(mn_mask_m))
      # for i in 1:length(mn_mask_m)
      #   println("starting $i")
      #   A = randn(mn_mask_m[i], n)
      #   x = randn(n)
      #   b = zeros(mn_mask_m[i])
      #   # Axtime_mul[i] = @elapsed mul!(b, A, x)
      #   Axtime_blas[i] = @elapsed BLAS.gemv!('N', true, A, x, false, b)
      #   println("finished $i")
      # end
      # axx.loglog(mn_mask_m, Axtime_blas, linestyle=":", label=L"$Ax$", color="red", zorder=1, linewidth=3)
      
      # legends
      # axx.get_yaxis().set_visible(false)
      # axx.add_artist(leg)
      # ax.get_yaxis().set_visible(true)
      ax.set_ylabel("time [s]")
      # ax.set_ylabel("scaling")
      println(PROJPATH * "figures/$(plotname)__$(n)__$(time_type).pdf")
      fig.tight_layout()
      fig.savefig(PROJPATH * "figures/$(plotname)__$(n)__$(time_type).pdf", bbox_inches = "tight", pad_inches=0.05)
    end
  end

  # fig2,ax2 = subplots()
  # ax2.cla()
  # Axtime_mul = zeros(5)
  # Axtime_blas = zeros(5)
  # for i in 1:5
  #   println("starting $i")
  #   A = randn(mn_mask_m[i], 100)
  #   x = randn(100)
  #   b = zeros(mn_mask_m[i])
  #   Axtime_mul[i] = @elapsed mul!(b, A, x)
  #   Axtime_blas[i] = @elapsed BLAS.gemv!('N', true, A, x, false, b)
  #   println("finished $i")
  # end
  # # ax22 = ax2.twinx()
  # ax22.cla()
  # ax22.loglog(mn_mask_m, mn_mask_m.^2, color="black", linewidth=3, linestyle=":")
  # ax22.loglog(mn_mask_m, mn_mask_m, color="black", linewidth=3, linestyle=":")
  # ax2.cla()
  # ax2.loglog(mn_mask_m, [sum(sum(out[method][mn_key]["time_A"])) for mn_key in mn_keys])
  # ax2.loglog(mn_mask_m, [out[method][mn_key]["time_A"][1][end] for mn_key in mn_keys])
  # ax2.loglog(mn_mask_m, [sum(sum(out[method][mn_key][f])) ./ sum(sum(out[method][mn_key]["lsit"])) for mn_key in mn_keys])
  # ax2.loglog(mn_mask_m, [out[method][mn_key]["time_A"][1][1] for mn_key in mn_keys])
  # ax2.loglog(mn_mask_m, Axtime_mul, color="red", linewidth=2)
  # ax2.loglog(mn_mask_m, Axtime_blas, color="blue", linewidth=2)
end

#=
fig2,ax2 = subplots()
ax2.cla()
axx2.cla()
axx2 = ax2.twinx()
ax2.plot(1:10, 2:11, label="main", zorder=100)
leg2 = ax2.legend(loc="center left", borderaxespad=1.)
leg2.remove()
axx2.plot(1:10, [100 for i in 1:10], label="second", zorder=10, color="black", linewidth=2)
# lines, labs = ax2.get_legend_handles_labels()
# lines2, labs2 = axx2.get_legend_handles_labels()
# axx2.legend(lines, labs, framealpha=1.0, loc="center left")
axx2.legend(loc="center left", borderaxespad=1.)
axx2.add_artist(leg2)
ax2.grid(true)
ax2.set_ylabel("main")
axx2.get_yaxis().set_visible(false)
fig2.savefig("/Users/jakeroth/Downloads/test.pdf")
close(fig2)


fig2, ax1 = plt.subplots()
fig.set_size_inches(18/1.5, 10/1.5)
ax2 = ax1.twinx()
ax1.plot(1:10, 1:10, label="main", linestyle="-")
ax1.grid(true)
ax2.plot(1:10, [100 for i in 1:10], label="second", linestyle="--")
# ax1.legend(loc="center left", borderaxespad=1.)
# ax2.legend(loc="upper right", borderaxespad=1.)
legend_1 = ax1.legend(loc="center left", borderaxespad=1.)
legend_1.remove()
# ax2.legend(loc="center left", borderaxespad=1.)
ax2.add_artist(legend_1)
fig2.savefig("/Users/jakeroth/Downloads/test.pdf")
close(fig2)
=#