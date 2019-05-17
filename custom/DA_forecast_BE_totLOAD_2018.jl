using RepresentativeDaysFinders
using JuMP
using PyCall
# using GLPK
using Gurobi
using CSV
using DataFrames
using Dates
using Timezones

##################################################################################
# collect data
##################################################################################
#done externally

##################################################################################
# preprocess data if needed
##################################################################################
#interpolation done here
df_LF_wind = ENTSOEcsv2dataframe(joinpath(@__DIR__, "windforecast_ELIA.csv"), 12, :LF_wind; delim = ';')
df_LF_wind_interpolated = interpolatedataframe(df_LF_wind, :LF_wind, 0.25,1.)
df_LF_solar = ENTSOEcsv2dataframe(joinpath(@__DIR__, "solarforecast_ELIA.csv"), 2, :LF_solar; delim = ';')
df_LF_solar_interpolated = interpolatedataframe(df_LF_solar, :LF_solar, 0.25,1.)
csv_wind_raw = CSV.read(joinpath(@__DIR__, "windforecast_ELIA.csv"), delim = ";")
csv_solar_raw = CSV.read(joinpath(@__DIR__, "solarforecast_ELIA.csv"), delim = ";")
df_test = join(csv_wind_raw, csv_solar_raw, on = [:Column1, :Column1], makeunique=true )


test1 = RepresentativeDaysFinders.check_temporal_consistency(csv_wind_raw, 0.25,"y-m-d HH:MM:SSz", :Column1)

df = DateFormat("y-m-d HH:MM:SSz")
# df_str = "y-m-d HH:MM:SSz"

# put all date in same time zone and add column to df
dt_df = Dates.Minute(15)
tzone = tz"Europe/Brussels"
csv_wind_raw[:date_formated] = map(str->astimezone(ZonedDateTime(str, df), tzone),csv_wind_raw[:Column1])
empty_row = similar(csv_wind_raw, 1)
for n in names(empty_row)
    if n != :date_formated
        empty_row[n] = missing
    end
end


for i in 2:length(csv_wind_raw[:date_formated])
    t_current = csv_wind_raw[:date_formated][i]
    t_prev = csv_wind_raw[:date_formated][i-1]
    dt =  Dates.Minute(t_current -  t_prev)
    if  dt != dt_df
        println((i, dates_date[i], dt))
        #add line
        closed = false
        t_add = t_prev
        while !closed
            tba = deepcopy(empty_row)
            t_add = t_add + dt_df
            tba[:date_formated] = t_add
            append!(csv_wind_raw, tba)
            if Dates.Minute(t_current - t_add) == dt_df
                closed = true
            end
        end
    end
end
sort(csv_wind_raw,:date_formated)
##################################################################################
# Specify location of config-file
##################################################################################
using Cbc
config_file = joinpath(@__DIR__,"DA_LOAD_forcast_BE_2018.yaml")
dft = findRepresentativeDays(config_file, with_optimizer(Cbc.Optimizer; seconds = 3600*5))
RepresentativeDaysFinders.create_plots(dft)
