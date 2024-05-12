function writeout_alm_grid(res_grid::Dict, datapath::String, PI::ProblemInstance, m::Vector, n::Integer, repi::Integer)
  ks = keys(res_grid)
  for k in ks
    if k != :walltime
      pct = res_grid[k][:k][1]/PI.m[1]
      println("pct = $pct")
      outname = "warm/k_$(Int(ceil(pct*1000)))pct__L_1"
      mkpath(datapath * "obj_$(PI.objtype)/$outname/")  
      writeout_alm(res_grid[k], datapath, PI, outname, m, n, repi)
    else
      writedlm(datapath * "obj_$(PI.objtype)/" * "grid_walltime.csv", res_grid[k])
    end
  end
end
function writeout_alm(res::Dict, datapath::String, PI::ProblemInstance, outname::String, m::Vector, n::Integer, repi::Integer)
  # common
  mkpath(datapath * "obj_$(PI.objtype)/$outname/")
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__pinfeas.csv", haskey(res, :pinfeas) ? res[:pinfeas] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__dinfeas.csv", haskey(res, :dinfeas) ? res[:dinfeas] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__relgap.csv", haskey(res, :relgap) ? res[:relgap] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__walltime.csv", haskey(res, :walltime) ? res[:walltime] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__walltime_i1.csv", haskey(res, :walltime_i1) ? res[:walltime_i1] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__walltime_i2.csv", haskey(res, :walltime_i2) ? res[:walltime_i2] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__pobj.csv", haskey(res, :pobj) ? res[:pobj] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__dobj.csv", haskey(res, :dobj) ? res[:dobj] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__iter.csv", haskey(res, :iter) ? res[:iter] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__status.csv", haskey(res, :status) ? res[:status] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__avg_step_size.csv", haskey(res, :avg_step_size) ? res[:avg_step_size] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__nnz.csv", sum(prod.(size.(PI.A))))
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__kktres.csv", haskey(res, :err) ? res[:err] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__r.csv", haskey(res, :maxksumrhs) ? res[:maxksumrhs] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__mks.csv", haskey(res, :maxksum) ? res[:maxksum] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__cvar.csv", haskey(res, :cvar) ? res[:cvar] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__cvar_infeas.csv", haskey(res, :cvar_infeas) ? res[:cvar_infeas] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__cvar_rhs.csv", haskey(res, :cvar_rhs) ? res[:cvar_rhs] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__normx.csv", haskey(res, :normx) ? res[:normx] : NaN)
  # specific
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__retcode.csv", haskey(res, :retcode) ? res[:retcode] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__ssnstatus.csv", haskey(res, :substatus) ? res[:substatus] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__feval.csv", haskey(res, :lsit) ? sum(res[:lsit]) : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__vphieval.csv", haskey(res, :lsit) ? sum(res[:lsit]) : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__gphieval.csv", haskey(res, :nit) ? sum(res[:nit]) : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__Hphieval.csv", haskey(res, :nit) ? sum(res[:nit]) : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__avg_nd_size.csv", haskey(res, :avg_nd_size) ? res[:avg_nd_size] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__ssnpinfeas.csv", haskey(res, :sub_pinfeas) ? res[:sub_pinfeas] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__ssndinfeas.csv", haskey(res, :sub_dinfeas) ? res[:sub_dinfeas] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__ssnrelgap.csv", haskey(res, :sub_relgap) ? res[:sub_relgap] : NaN)
  # timing
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_proj.csv", haskey(res, :ctime_proj) ? res[:ctime_proj] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_sort.csv", haskey(res, :ctime_sort) ? res[:ctime_sort] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_nd.csv", haskey(res, :ctime_nd) ? res[:ctime_nd] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_Hphi.csv", haskey(res, :ctime_Hphi) ? res[:ctime_Hphi] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_gphi.csv", haskey(res, :ctime_gphi) ? res[:ctime_gphi] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_vphi.csv", haskey(res, :ctime_vphi) ? res[:ctime_vphi] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_dirn.csv", (haskey(res, :ctime_nd) && haskey(res, :ctime_Hphi)) ? res[:ctime_nd] + res[:ctime_Hphi] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_A.csv", haskey(res, :ctime_Ax) ? res[:ctime_Ax] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_At.csv", haskey(res, :ctime_Atlam) ? res[:ctime_Atlam] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__time_ssn.csv", haskey(res, :timessn) ? res[:timessn] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__lsit.csv", haskey(res, :lsit) ? res[:lsit] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__sps.csv", haskey(res, :sps) ? res[:sps] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__psc.csv", haskey(res, :partial_sort_count) ? res[:partial_sort_count] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__nit.csv", haskey(res, :nit) ? res[:nit] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__sigma.csv", haskey(res, :sigma) ? res[:sigma] : NaN)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__alm__tau.csv", haskey(res, :tau) ? res[:tau] : NaN)
end
function writeout_grb(res::Dict, datapath::String, PI::ProblemInstance, outname::String, m::Vector, n::Integer, repi::Integer)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__pinfeas.csv", res[:pinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__dinfeas.csv", res[:dinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__relgap.csv", res[:relgap])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__walltime.csv", res[:walltime])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__pobj.csv", res[:pobj])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__iter.csv", res[:iter])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__status.csv", res[:status])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__t.csv", res[:t])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__nnz.csv", res[:nnz])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__mks.csv", res[:maxksum])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__r.csv", res[:maxksumrhs])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__cvar.csv", res[:cvar]) # vector
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__cvar_infeas.csv", res[:cvar_infeas]) # scalar
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__cvar_rhs.csv", res[:cvar_rhs]) # vector
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb__normx.csv", res[:normx]) # scalar
end
function writeout_grb_oa(res::Dict, datapath::String, PI::ProblemInstance, outname::String, m::Vector, n::Integer, repi::Integer)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__pinfeas.csv", res[:pinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__dinfeas.csv", res[:dinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__relgap.csv", res[:relgap])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__walltime.csv", res[:walltime])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__subtimes.csv", res[:times])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__subiters.csv", res[:iters])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__pobj.csv", res[:pobj])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__iter.csv", res[:iter])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__status.csv", res[:status])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__ncon.csv", res[:numcon])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__mks.csv", res[:maxksum])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__r.csv", res[:maxksumrhs])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__cvar.csv", res[:cvar]) # vector
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__cvar_infeas.csv", res[:cvar_infeas]) # scalar
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__cvar_rhs.csv", res[:cvar_rhs]) # vector
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_oa__normx.csv", res[:normx]) # scalar
end
function writeout_grb_qr(res::Dict, datapath::String, PI::ProblemInstance, outname::String, m::Vector, n::Integer, repi::Integer)
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_qr__pinfeas.csv", res[:pinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_qr__dinfeas.csv", res[:dinfeas])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_qr__relgap.csv", res[:relgap])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_qr__walltime.csv", res[:walltime])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_qr__pobj.csv", res[:pobj])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_qr__iter.csv", res[:iter])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_qr__status.csv", res[:status])
  writedlm(datapath * "obj_$(PI.objtype)/$outname/$(m[1])_$(n)_$(repi)__grb_qr__normx.csv", res[:normx]) # scalar
end

#
# get regression data
#

function get_regression_data(fpath::String, fname::String, loglog::Bool=true, stdize::Bool=true, writeout::Bool=false)
  ds = Parquet2.Dataset(fpath * fname) # load metadata
  Parquet2.appendall!(ds) # load metadata
  df = DataFrame(ds; copycols=false) # read data (~40GB)

  # categorical columns
  catcols = [
    "UBER",
    "PUborough_EWR",
    "PUborough_QNS",
    "PUborough_BNX",
    "PUborough_MTN",
    "PUborough_STI",
    "PUborough_BKN",
    "PUborough_UNK",
    "DOborough_EWR",
    "DOborough_QNS",
    "DOborough_BNX",
    "DOborough_MTN",
    "DOborough_STI",
    "DOborough_BKN",
    "DOborough_UNK",
    "shared_match_flag_Y",
    "shared_request_flag_Y",
    "wav_match_flag_Y",
    "wav_request_flag_Y",
    "tod_early",
    "tod_morning",
    "tod_afternoon",
    "tod_evening",
    "tod_night",
    "tod_late",
    "light_out",
    "weekend",
    "cloudy",
    "congestion_Y",
  ]

  # drop helper columns
  dropcols = [
    # helper
    "hvfhs_license_num",
    "dispatching_base_num",
    "originating_base_num",
    "request_datetime",
    "on_scene_datetime",
    "pickup_datetime",
    "dropoff_datetime",
    "PULocationID",
    "DOLocationID",
    "shared_request_flag",
    "shared_match_flag",
    "access_a_ride_flag",
    "wav_request_flag",
    "wav_match_flag",
    "PUborough",
    "DOborough",
    "JUNO",
    "VIA",
    "datetime_unix",
    "weather_row",
    "PUborough_MTN", # default PU borough = Manhattan
    "DOborough_MTN", # default DO borough = Manhattan
    "tod_afternoon", # default tod = afternoon
    "LYFT", # default carrier = LYFT
    "clouds",
    # nonhelpful
    "airport_fee",
    "bcf",
    "congestion_surcharge",
    "sales_tax",
    "snow1h",
    "snow3h",
    "wind_deg",
    "light_out",
  ]
  select!(df, Not(dropcols))

  # loglog
  if loglog
    for n in names(df)
      println(n, "  ", typeof(df[!,n]))
      if isa(df[!,n], BitVector) || n ∈ catcols
        continue
      end
      if any(df[!,n] .== -1.0)
        println("!ERR!")
      end
      df[!,n] .= log.(1.0 .+ df[!,n])
    end
  end
  GC.gc()
  
  # standardize
  if stdize
    for n in names(df)
      println(n, "  ", typeof(df[!,n]))
      if isa(df[!,n], BitVector) || n ∈ catcols
        continue
      end
      df[!,n] .= (df[!,n] .- mean(df[!,n])) ./ std(df[!,n])
    end
  end
  GC.gc()

  # populate matrix
  dependent_name = "tips"
  Amat = zeros(Float64, size(df,1), size(df,2)-1)
  yvec = zeros(Float64, size(df,1))
  yvec .= df.tips
  select!(df, Not(dependent_name))
  for (i,n) in enumerate(names(df))
    Amat[:,i] .= df[!,n]
  end

  if writeout
    fname_name, fname_type = split(fname, ".")
    df[:, dependent_name] = yvec
    Parquet2.writefile(fpath * fname_name * "_standardized." * fname_type, df)
  end
  GC.gc()
  return (Amat, yvec, names(df))
end