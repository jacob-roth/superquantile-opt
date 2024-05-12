import Pkg; Pkg.activate("."); Pkg.instantiate()
using DelimitedFiles, Parquet2, DataFrames
using HTTP, JSON # Downloads
using StatsBase, Dates
global const DATAURL = "https://d37ci6vzurychx.cloudfront.net/trip-data/"
global const DATAPATH = "/home/roth0674/drive/taxidata/"
global const DATATYPE = "fhvhv"
global const NYCLL = (40.7826, -73.9656)
global const APIKEY = "TO_BE_PROVIDED"

#
# set time period
#

const year_0 = 2022
const month_0 = 10
const nmonths = 12
const nhours = 365 * 24
months = [DateTime("$(year_0)-$(month_0)", "yyyy-mm") + i*Month(1) for i in 0:nmonths-1]
hours = [DateTime("$(year_0)-$(month_0)-00", "yyyy-mm-HH") + i*Hour(1) for i in 0:nhours-1]

#
# get taxi data
#

function get_data_month(month::Integer, year::Integer, datatype::String=DATATYPE, dataurl::String=DATAURL)
  """
  datatype ∈ {"fhvhv" := high volume for-hire-vehicle, "fhv" := for-hire-vehicle, "green" := green taxi, "yellow" := yellow taxi}
  note: "fhvhv" has the most data fields, dates back to 2019-02
  """
  fname = "$(datatype)_tripdata_$(year)-$(lpad(month, 2, "0")).parquet"
  loc = datapath * "parquet/$fname"
  url = dataurl * fname
  myfile = Downloads.download(url, loc)
end

for dt in months
  println(dt)
  get_data_month(Month(dt).value, Year(dt).value)
end

#
# get weather data
#

function apistring(lat::Float64, lon::Float64, dt::DateTime, apikey::String)
  dt_int = round(Int, datetime2unix(dt))
  return "https://api.openweathermap.org/data/3.0/onecall/timemachine?lat=$(lat)&lon=$(lon)&dt=$(dt_int)&appid=$(apikey)"
end

for (i,dth) in enumerate(hourds)
  println("$dth | $i / $(length(hours)) = " * @sprintf("%0.4f", i / length(hours)))
  r = HTTP.get(apistring(NYCLL[1], NYCLL[2], dth, APIKEY))
  d = JSON.parse(String(r.body))
  open(DATAPATH * "weather/$(round(Int, datetime2unix(dth))).json","w") do f
    JSON.print(f, d)
  end
end

weatherdir = DATAPATH * "weather/"
weatherfields = [
  "dt",          #  1
  "sunrise",     #  2
  "sunset",      #  3
  # -----------------
  "visibility",  #  4
  "dew_point",   #  5
  "wind_deg",    #  6
  "wind_speed",  #  7
  "wind_gust",   #  8
  "temp",        #  9
  "feels_like",  # 10
  "humidity",    # 11
  "pressure",    # 12
  "clouds",      # 13
  # -----------------
  "rain1h",      # 14
  "rain3h",      # 15
  "snow1h",      # 16
  "snow3h",      # 17
  # -----------------
  "light_out", # vs dark
  "weekend", # vs weekday
]
weatherfieldsdict = Dict()
for f in weatherfields
  weatherfieldsdict[f] = findfirst(f .== weatherfields)
end
jsons = readdir(weatherdir)[[occursin(".json", x) for x in readdir(weatherdir)]]
weatherdata = zeros(Float64, length(jsons), length(weatherfields))
for (i,f) in enumerate(jsons)
  d = JSON.parsefile(weatherdir * f)["data"][1]
  for j in eachindex(weatherfields)
    if j <= 13
      if haskey(d, weatherfields[j])
        weatherdata[i,j] = d[weatherfields[j]]
      end
    elseif j == 14 # rain1h
      if haskey(d, "rain")
        weatherdata[i,j] = haskey(d["rain"], "1h") ? d["rain"]["1h"] : 0.0
      end
    elseif j == 15 # rain3h
      if haskey(d, "rain")
        weatherdata[i,j] = haskey(d["rain"], "3h") ? d["rain"]["3h"] : 0.0
      end
    elseif j == 16 # snow1h
      if haskey(d, "snow")
        weatherdata[i,j] = haskey(d["snow"], "1h") ? d["snow"]["1h"] : 0.0
      end
    elseif j == 17 # snow1h
      if haskey(d, "snow")
        weatherdata[i,j] = haskey(d["snow"], "3h") ? d["snow"]["3h"] : 0.0
      end
    elseif j == 18 # light_out
      weatherdata[i,j] = (weatherdata[i,1] >= weatherdata[i,2]) && (weatherdata[i,1] <= weatherdata[i,3])
    elseif j == 19 # weekend
      weatherdata[i,j] = (dayname(unix2datetime(weatherdata[i,1])) == "Saturday") || (dayname(unix2datetime(weatherdata[i,1])) == "Sunday")
    end
  end
end
writedlm(weatherdir * "weatherdata.csv", weatherdata)
writedlm(weatherdir * "weatherfields.csv", weatherfields)

#
# load taxidata
#

# load taxidata
taxidir = DATAPATH * "parquet/"
tds = Parquet2.Dataset(taxidir) # load metadata
Parquet2.appendall!(tds) # load metadata
tdf = DataFrame(tds; copycols=false) # read data (40GB)
mask = tdf[:,"driver_pay"] .> 0 # pay>0
mask .*= tdf[:,"tips"] .> 0 # used CC to tip
mask .*= tdf[:,"base_passenger_fare"] .> 0 # fare>0
mask .*= (tdf[:,"trip_miles"] .> 0) .* (tdf[:,"trip_miles"] .<= 50) # 0<miles<=50mi
mask .*= (tdf[:,"trip_time"] .> 0) .* (tdf[:,"trip_time"] .<= 6*60*60) # 0<time<=6h
tdf = tdf[mask, All()] # 55GB
GC.gc()
taxi_zone_lookup = readdlm(DATAPATH * "taxi_zone_lookup.csv", ',')

# missing
mask = .!(ismissing.(tdf.on_scene_datetime) .|| ismissing.(tdf.request_datetime) .|| ismissing.(tdf.pickup_datetime)) # nonmissing
tdf = tdf[mask, All()]
GC.gc()

# physical timing
mask = (tdf.request_datetime .<= tdf.on_scene_datetime) .&& (tdf.on_scene_datetime .<= tdf.pickup_datetime) .&& (tdf.pickup_datetime .<= tdf.dropoff_datetime) # timing
tdf = tdf[mask, All()]
GC.gc()

# compute: delay, speed
tdf[:, "delay_1"] = tdf[!, "on_scene_datetime"] .- tdf[!, "request_datetime"] # milliseconds
tdf[:, "delay_2"] = tdf[!, "pickup_datetime"] .- tdf[!, "on_scene_datetime"] # milliseconds
tdf[!, "delay_1"] .= [x.value/(1000*60) for x in tdf[!, "delay_1"]] # min
tdf[!, "delay_2"] .= [x.value/(1000*60) for x in tdf[!, "delay_2"]] # min
# tdf[!, "delay_1_squared"] .= tdf[!, "delay_1"].^2 #! square after standardizing
# tdf[!, "delay_2_squared"] .= tdf[!, "delay_2"].^2 #! square after standardizing
tdf[:, "avg_speed"] = tdf[!, "trip_miles"] ./ tdf[!, "trip_time"] * 60 * 60 # mph
tdf = tdf[findall(tdf[:, "avg_speed"] .<= 70), All()]
GC.gc()

# fare per
tdf[:, "fare_per_mile"] = tdf[!, "base_passenger_fare"] ./ tdf[!, "trip_miles"]
tdf[:, "fare_per_min"] = tdf[!, "base_passenger_fare"] ./ (tdf[!, "trip_time"]./60)
tdf[:, "pay_per_mile"] = tdf[!, "driver_pay"] ./ tdf[!, "trip_miles"]
tdf[:, "pay_per_min"] = tdf[!, "driver_pay"] ./ (tdf[!, "trip_time"]./60)
tdf[:, "driver_share"] = tdf[!, "driver_pay"] ./ tdf[!, "base_passenger_fare"]
tdf[:, "congestion_per_min"] = tdf[!, "congestion_surcharge"] ./ (tdf[!, "trip_time"]./60)
tdf[:, "congestion_per_mile"] = tdf[!, "congestion_surcharge"] ./ (tdf[!, "trip_miles"])

# one-hot
tdf[:, "JUNO"] = [x .== "HV0002" for x in tdf[:, "hvfhs_license_num"]] # none
tdf[:, "UBER"] = [x .== "HV0003" for x in tdf[:, "hvfhs_license_num"]]
tdf[:, "VIA"] = [x .== "HV0004" for x in tdf[:, "hvfhs_license_num"]] # none
tdf[:, "LYFT"] = [x .== "HV0005" for x in tdf[:, "hvfhs_license_num"]]
tdf[:, "PUborough"] = [taxi_zone_lookup[ID, 2] for ID in tdf[:, "PULocationID"]]
tdf[:, "DOborough"] = [taxi_zone_lookup[ID, 2] for ID in tdf[:, "DOLocationID"]]
tdf[:, "PUborough_EWR"] = tdf[:, "PUborough"] .== "EWR"
tdf[:, "PUborough_QNS"] = tdf[:, "PUborough"] .== "Queens"
tdf[:, "PUborough_BNX"] = tdf[:, "PUborough"] .== "Bronx"
tdf[:, "PUborough_MTN"] = tdf[:, "PUborough"] .== "Manhattan"
tdf[:, "PUborough_STI"] = tdf[:, "PUborough"] .== "Staten Island"
tdf[:, "PUborough_BKN"] = tdf[:, "PUborough"] .== "Brooklyn"
tdf[:, "PUborough_UNK"] = tdf[:, "PUborough"] .== "Unknown"
tdf[:, "DOborough"] = [taxi_zone_lookup[ID, 2] for ID in tdf[:, "DOLocationID"]]
tdf[:, "DOborough_EWR"] = tdf[:, "DOborough"] .== "EWR"
tdf[:, "DOborough_QNS"] = tdf[:, "DOborough"] .== "Queens"
tdf[:, "DOborough_BNX"] = tdf[:, "DOborough"] .== "Bronx"
tdf[:, "DOborough_MTN"] = tdf[:, "DOborough"] .== "Manhattan"
tdf[:, "DOborough_STI"] = tdf[:, "DOborough"] .== "Staten Island"
tdf[:, "DOborough_BKN"] = tdf[:, "DOborough"] .== "Brooklyn"
tdf[:, "DOborough_UNK"] = tdf[:, "DOborough"] .== "Unknown"
tdf[:, "shared_match_flag_Y"] = [x == "Y" for x in tdf[:, "shared_match_flag"]]
tdf[:, "shared_request_flag_Y"] = [x == "Y" for x in tdf[:, "shared_request_flag"]]
tdf[:, "wav_match_flag_Y"] = [x == "Y" for x in tdf[:, "wav_match_flag"]]
tdf[:, "wav_request_flag_Y"] = [x == "Y" for x in tdf[:, "wav_request_flag"]]
tdf[!, "tod_early"] = (hour.(tdf.on_scene_datetime) .>= 4) .&& (hour.(tdf.on_scene_datetime) .< 7) # 4 - 7am
tdf[!, "tod_morning"] = (hour.(tdf.on_scene_datetime) .>= 7) .&& (hour.(tdf.on_scene_datetime) .< 11) # 7 - 11am
tdf[!, "tod_afternoon"] = (hour.(tdf.on_scene_datetime) .>= 11) .&& (hour.(tdf.on_scene_datetime) .< 15) # 11am - 3pm
tdf[!, "tod_evening"] = (hour.(tdf.on_scene_datetime) .>= 15) .&& (hour.(tdf.on_scene_datetime) .< 19) # 3pm - 7pm
tdf[!, "tod_night"] = (hour.(tdf.on_scene_datetime) .>= 19) .&& (hour.(tdf.on_scene_datetime) .< 23) # 7pm - 11pm
tdf[!, "tod_late"] = (hour.(tdf.on_scene_datetime) .>= 23) .|| (hour.(tdf.on_scene_datetime) .< 4) # 11pm - 4am
GC.gc()

#
# merge datasets
#

# get weatherdata
weatherdir = DATAPATH * "weather/"
weatherdata = readdlm(weatherdir * "weatherdata.csv")

# add weatherdata (5-10min)
dlo = unix2datetime.(weatherdata[1:end,1])
dhi = [unix2datetime.(weatherdata[2:end,1]); unix2datetime(weatherdata[end,1]) + Hour(1)]
tdf[:, "weather_row"] = [findfirst((dt .> dlo) .* (dt .<= dhi)) for dt in tdf[:, "request_datetime"]]
tdf[!, "weather_row"][isnothing.(tdf[:, "weather_row"])] .= 1 #! assign first weather to rides that finish in the time window but start before it
tdf[!, "weather_row"] = convert(Vector{Int64}, tdf[!, "weather_row"])
tdf[:, "visibility"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["visibility"]] for i in 1:nrow(tdf)]
tdf[:, "dew_point"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["dew_point"]] for i in 1:nrow(tdf)]
tdf[:, "wind_deg"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["wind_deg"]] for i in 1:nrow(tdf)]
tdf[:, "wind_speed"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["wind_speed"]] for i in 1:nrow(tdf)]
tdf[:, "wind_gust"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["wind_gust"]] for i in 1:nrow(tdf)]
tdf[:, "temp"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["temp"]] for i in 1:nrow(tdf)]
tdf[:, "feels_like"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["feels_like"]] for i in 1:nrow(tdf)]
tdf[:, "humidity"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["humidity"]] for i in 1:nrow(tdf)]
tdf[:, "pressure"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["pressure"]] for i in 1:nrow(tdf)]
tdf[:, "clouds"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["clouds"]] for i in 1:nrow(tdf)]
tdf[:, "rain1h"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["rain1h"]] for i in 1:nrow(tdf)]
tdf[:, "rain3h"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["rain3h"]] for i in 1:nrow(tdf)]
tdf[:, "snow1h"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["snow1h"]] for i in 1:nrow(tdf)]
tdf[:, "snow3h"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["snow3h"]] for i in 1:nrow(tdf)]
tdf[:, "light_out"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["light_out"]] for i in 1:nrow(tdf)]
tdf[:, "weekend"] = [weatherdata[tdf[i, "weather_row"], weatherfieldsdict["weekend"]] for i in 1:nrow(tdf)]
tdf[:, "cloudy"] =  tdf[:, "clouds"] .> 50
tdf[:, "congestion_Y"] =  tdf[:, "congestion_surcharge"] .> 0
GC.gc()
#! 80GB

fpath = DATAPATH
fname = DATAFILENAME
ds = Parquet2.Dataset(fpath * fname) # load metadata
Parquet2.appendall!(ds) # load metadata
tdf = DataFrame(ds; copycols=false) # read data (~40GB)
tdf[:,"congestion_surcharge_pwl"] = (tdf[:,"congestion_surcharge"] .- minimum(tdf[:,"congestion_surcharge"])) .* (tdf[:,"congestion_surcharge"] .> minimum(tdf[:,"congestion_surcharge"]))


#
# further sanity checking
#

[quantile(tdf[!,"fare_per_mile"], 0.999+0.0001i) for i in 1:9]
[quantile(tdf[!,"fare_per_min"], 0.9999999+0.00000001i) for i in 1:9]
mask = tdf[!,"fare_per_min"] .<= 30 # half of https://www.businessinsider.com/ubers-highest-surge-price-ever-may-be-50x-2014-11?utm_source=copy-link&utm_medium=referral&utm_content=topbar
mask .*= .!((tdf[!,"fare_per_mile"] .> 100) .&& (tdf[!,"trip_time"] .> (30*60)))
mask .*= tdf[!,"pay_per_min"] .<= 30 # half of https://www.businessinsider.com/ubers-highest-surge-price-ever-may-be-50x-2014-11?utm_source=copy-link&utm_medium=referral&utm_content=topbar
mask .*= .!((tdf[!,"pay_per_mile"] .> 100) .&& (tdf[!,"trip_time"] .> (30*60)))
tdf = tdf[mask,All()]
GC.gc()
mask = tdf[!, "driver_pay"].>1
[quantile(tdf[!,"driver_share"], 0.999+0.0001i) for i in 1:9]

#
# plot single variable relationships
#

# idx = 1:nrow(tdf)
# ridx=rand(idx, 5000)
for n in names(tdf)
  println(n)
  if isfile("./scatters/$n.pdf")
    continue
  else
    ax.cla();ax.scatter(
      (tdf[ridx, n]),
      (tdf[ridx, "tips"]), alpha=0.1
    );fig.savefig("./scatters/$n.pdf")
  end
end
for n in names(tdf)
  println(n)
  if false#isfile("./loglogscatters/$n.pdf")
    continue
  else
    if eltype(tdf[ridx, n])==String||eltype(tdf[ridx, n])==DateTime||occursin("log",n)
      continue
    else
      ax.cla();ax.scatter(
        log.(1.0.+tdf[ridx, n]),
        log.(tdf[ridx, "tips"]), alpha=0.1
      );fig.savefig("./loglogscatters/$n.pdf")
    end
  end
end

#
# interactions
#

tdf[!,"congestion_Y"] .= tdf[!,"congestion_surcharge"] .> 0 # https://stats.stackexchange.com/a/184371/113279

#
# check colinearity
#

ax.cla();ax.scatter(
  (tdf[ridx, "base_passenger_fare"]),
  (tdf[ridx, "bcf"]), alpha=0.1
);fig.savefig("./loglogscatters/z_colinearity.pdf") # little variablility; 

#=
SUMMARY

====== ~linear on linear-linear ======
- base_passenger_fare
- bcf
- driver_pay
- rain1h (perhaps; perhaps exponential relationship)
- rain3h (perhaps; perhaps exponential relationship)
- sales_tax
- trip_miles
- trip_time

====== not linear on linear-linear ======
- avg_speed
- congestion_per_mile (roughly exponential)
- congestion_per_min (roughly exponential)
- delay_1
- delay_2
- dew_point
- driver_share
- fare_per_mile (roughly exponential)
- fare_per_min
- feels_like
- humidity
- pay_per_mile (roughly exponential)
- pressure
- temp
- temp_squared_centered (roughly exponential?)
- tolls
- visibility
- wind_gust
- wind_speed
====== remove linear-linear ======
- pay_per_min
- snow1h (not enough nonzero)
- snow3h (not enough nonzero)
- wind_deg

====== ~linear on log-log ======
- avg_speed
- base_passenger_fare
- congestion_per_mile
- congestion_per_min
- driver_pay
- fare_per_mile
- fare_per_min
- pay_per_mile
- pay_per_min
- rain1h (maybe)
- rain3h (maybe)
- tolls (maybe)
- trip_miles
- trip_time
- wind_gust (maybe)
- visibility (maybe)
====== not linear on log-log ======
- delay_1 (none apparent)
- delay_2 (none apparent)
- dew_point (none)
- driver_share (none)
- humidity (none)
- feels_like (none)
- pressure (none)
- temp (none)
- wind_speed (none)
====== remove log-log ======
- airport_fee (roughly discrete; handle by borough)
- bcf (proportion of fare? different proportions over time?)
- congestion_surcharge (roughly discrete)
- sales_tax (artifact at zero tax)
- snow1h (not enough nonzero)
- snow3h (not enough nonzero)
- wind_deg
- clouds

====== roughly colinear ======
- sales_tax and base_passenger_fare
- bcf and base_passenger_fare
=#

#
# write out
#

combineddir = DATAPATH * "parquet_combined2/"
for n in names(tdf);
  println(n)
  if eltype(first(tdf[!,n])) == Char
    tdf[!,n] .= convert(Vector{String}, tdf[!,n])
  else
    tdf[!,n] .= convert(Vector{eltype(first(tdf[!,n]))}, tdf[!,n])
  end
end
# Parquet2.writefile(combineddir * "taxi_weather_data_nomissing.parquet", tdf)
Parquet2.writefile("/home/roth0674/drive/taxidata/parquet_combined2/taxi_weather_data_nomissing.parquet", tdf)
