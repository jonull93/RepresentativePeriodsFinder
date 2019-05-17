using PyCall
datetime = pyimport("datetime")
pd = pyimport("pandas")
sc = pyimport("elia_xml_services.sources")
Region = pyimport("elia_xml_services.helpers").Region
WindConnection = pyimport("elia_xml_services.helpers").WindConnection
WindLocation = pyimport("elia_xml_services.helpers").WindLocation

year = 2017
start_period = datetime.datetime(2017, 1, 1)
end_period = datetime.datetime(2017, 12, 31)

wf = sc.WindForecast()
sf = sc.SolarForecast()

sdf = sf.data(start_period, end_period, region=Region.Limburg)
wdf = wf.data(start_period, end_period, location=WindLocation.Both, connection=WindConnection.Elia_Connected)

# write data
wdf.to_csv(joinpath(@__DIR__, "windforecast.csv"), sep=";")
sdf.to_csv(joinpath(@__DIR__, "solarforecast.csv"), sep=";")
