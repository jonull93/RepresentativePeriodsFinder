using RepresentativePeriodsFinder
using YAML
using CSV
using Gurobi, JuMP, Cbc
# using CSVFiles
using DataFrames
using XLSX
RPF = RepresentativePeriodsFinder
# using Dates

function repr_year_data(inp_name)

    # inp_name = "2018"
    # Create PeriodsFinder
    cf_name="periods_15_19_$(inp_name).yaml"
    config_file = normpath(joinpath(@__DIR__, "periods_15_19_undef.yaml"))
    cf = YAML.load_file(config_file)

    cf["time_series"]["p$(inp_name)"] = Dict()
    d_path= "input_data/$(inp_name).csv"
    cf["time_series"]["default"]["source"] = d_path
    cf["time_series"]["p$(inp_name)"]["source"] = d_path
    cf["time_series"]["p$(inp_name)"]["csv_options"] = Dict()
    cf["time_series"]["p$(inp_name)"]["csv_options"]["delim"] = "," 
    # cf["time_series"]["p$(inp_name)"]["csv_options"]["dateformat"] = "" 
    # cf["time_series"]["p$(inp_name)"]["csv_options"]["dateformat"] = "yyyy-mm-dd HH:MM:SS" 
    cf["time_series"]["p$(inp_name)"]["timestamp"] = "date" 
    cf["time_series"]["p$(inp_name)"]["timestamp_column"] = "date" 
    cf["time_series"]["p$(inp_name)"]["value_column"] = "load" 
    cf["time_series"]["p$(inp_name)"]["weight"] = 1 
    # cf["time_series"]["p$(inp_name)"]["interpolation_type"] = "" 
    cf["time_series"]["p$(inp_name)"]["sampling_time"] = "Hour(1)" 

    cf["results"]["result_dir"] = "results_$(inp_name)"

    # reformat_date(d_path)

    dpro=DataFrame(CSV.File(joinpath(@__DIR__,d_path)))
    imax= argmax(dpro[:,:load])

    cf["method"]["options"]["mandatory_periods"] = [imax]

    cf_path = normpath(joinpath(@__DIR__, cf_name))
    YAML.write_file(cf_path, cf)


    pf = PeriodsFinder(cf_path, populate_entries=true)

    # Delete entries for time series and duration curve error
    delete!(pf.config["method"]["optimization"], "time_series_error")
    delete!(pf.config["method"]["options"], "ordering_error")
    delete!(pf.config["method"], "clustering")

    # Make life a bit easier for the optimiser, reduce the number of days
    N_total = 8759
    N_rep = 20
    pf.config["method"]["options"]["total_periods"] = N_total
    pf.config["method"]["options"]["representative_periods"] = N_rep

    # Setup optimizer
    opt = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 30)

    pf.config["method"]["optimization"]["duration_curve_error"]["type"] = "absolute"
    find_representative_periods(pf, optimizer=opt, reset=true)
end

repr_year_data("2019")


# function reformat_date(p)
#     df = DataFrame(CSV.File(p))
#     datestring = df[:,:date]
#     datedate = [DateTime(i, DateFormat("yyyy-mm-dd HH:MM:SS")) for i in datestring]
#     df_new = DataFrame([datedate,df[:,:load]], [:date, :load])
#     CSV.write(p, df_new)
# end
# reformat_date("input_data/2019.csv")


save_name = "5_dems_dyn_weight"
inp_names= ["2015","2016","2017","2018","2019"]

#find repr periods
[repr_year_data(inp) for inp in inp_names]

out_names = ["results_$(i)" for i in inp_names]

res_all = []
for (i,o) in zip(inp_names,out_names) 
    ds= DataFrame(CSV.File(joinpath(@__DIR__, o, "decision_variables_short.csv")))[:,:weights]
    re= DataFrame(CSV.File(joinpath(@__DIR__, o, "resulting_profiles.csv")))[:,Symbol("p$i")]
    push!(res_all, (ds, re))
end

N_timesteps = 1
N_periods   = length(res_all[1][1])
N_years     = 1

tot_length = N_timesteps * N_periods * N_years

StartDateTime=["" for t in 1:tot_length]
Year        = ones(tot_length,1)
Period      = ones(N_periods,1).*[i for i in 1:N_periods]
TimeStep    = ones(N_periods,1).*[1 for i in 1:N_periods]

scenario_names = ["s$(lpad(s,2,"0"))" for s=1:length(out_names)]

ED_title     = ["dHourlyEnergyMarket;r01_z01;$(si)" for si in scenario_names]
weight_title = ["Weights;$(si)" for si in scenario_names]


function ci(i,j)
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    "$(uppercase(alphabet[j]))$i"
end

XLSX.openxlsx(joinpath(@__DIR__,"$(save_name).xlsx"), mode="w")  do xf
# xf = XLSX.openxlsx("eldest_frame.xlsx", mode="w")
    sheet = xf[1]
    XLSX.rename!(sheet, "TimeSteps")
    sheet["A1"] = "Year"
    sheet["A2"] = Year
    sheet["B1"] = "Period"
    sheet["B2"] = Period
    sheet["C1"] = "TimeStep"
    sheet["C2"] = TimeStep
    sheet[1,4] = "StartDateTime"
    # sheet[2,4] = ones(N_periods,1).*0
    sheet[1,5] = ED_title
    sheet[ci(1,5+length(ED_title))] = weight_title
    for (i,(ds,re)) in enumerate(res_all)
        @show i
        sheet[ci(2,5-1+i)] = reshape(re, N_periods,1)
        sheet[ci(2,5+length(ED_title)-1+i)] = reshape(ds, N_periods,1)
    end
end

