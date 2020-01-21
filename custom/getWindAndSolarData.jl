using PyCall
datetime = pyimport("datetime")
pd = pyimport("pandas")
sc = pyimport("elia_xml_services.sources")
Region = pyimport("elia_xml_services.helpers").Region
WindConnection = pyimport("elia_xml_services.helpers").WindConnection
WindLocation = pyimport("elia_xml_services.helpers").WindLocation

year = 2017
start_period = datetime.datetime(2017, 1, 1)
end_period = datetime.datetime(2018, 1, 1)

wf = sc.WindForecast()
sf = sc.SolarForecast()

sdf = sf.data(start_period, end_period, region=Region.Limburg)
wind_Offshore_df = wf.data(start_period, end_period, location=WindLocation.Offshore, connection=WindConnection.Elia_Connected)
wind_Onshore_df = wf.data(start_period, end_period, location=WindLocation.Onshore, connection=WindConnection.Elia_Connected)

# write data
sdf.to_csv(joinpath(@__DIR__, "solarforecast_new.csv"), sep=";")
wind_Offshore_df.to_csv(joinpath(@__DIR__, "windforecast_Offshore_new.csv"), sep=";")
wind_Onshore_df.to_csv(joinpath(@__DIR__, "windforecast_Onshore_new.csv"), sep=";")

# function try_and_solve()
#     # @show i = 3
#     # runit = true
#     # while runit
#     #     try
#     #     start_period = datetime.datetime(2017, i, 1)
#     #     end_period = datetime.datetime(2017, i+2, 1)
#     #     wind_Onshore_df = wf.data(start_period, end_period, location=WindLocation.Onshore, connection=WindConnection.Elia_Connected)
#     #     wind_Onshore_df.to_csv(joinpath(@__DIR__, "windforecast_Onshore_new_month$(i).csv"), sep=";")
#     #     i+=2
#     #     catch
#     #     end
#     #     if i >= 12
#     #         @show runit = false
#     #     end
#     # end
#     runit = true
#     while runit
#         start_period = datetime.datetime(2017, 11, 1)
#         end_period = datetime.datetime(2018, 1, 1)
#         wind_Onshore_df = wf.data(start_period, end_period, location=WindLocation.Onshore, connection=WindConnection.Elia_Connected)
#         wind_Onshore_df.to_csv(joinpath(@__DIR__, "windforecast_Onshore_new_month_last.csv"), sep=";")
#         runit = false
#     end
# end
# try_and_solve()
