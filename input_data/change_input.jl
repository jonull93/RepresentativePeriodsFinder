using JuMP
using ExcelReaders
using DataFrames
using CSV

cd("C:/Users/GOVAERTN/Desktop/RepresentativeDaysFinder/input_data")
data = CSV.read("Input.csv")

DEM = zeros(8760)
PV = zeros(8760)
WIND = zeros(8760)
for j = 1:8760
    DEM[j] = get(data[:DEM][4*(j-1)+1] + data[:DEM][4*(j-1)+2]
            	    + data[:DEM][4*(j-1)+3] + data[:DEM][4*(j-1)+4])
    PV[j] = get(data[:PV][4*(j-1)+1] + data[:PV][4*(j-1)+2]
                        + data[:PV][4*(j-1)+3] + data[:PV][4*(j-1)+4])/4
    WIND[j] = get(data[:WIND][4*(j-1)+1] + data[:WIND][4*(j-1)+2]
                        + data[:WIND][4*(j-1)+3] + data[:WIND][4*(j-1)+4])/4
end

d = DataFrame(DEM = DEM)
d[:PV] = PV
d[:WIND] = WIND

CSV.write("in.csv",d)
