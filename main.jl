using RepresentativeDaysFinder

##################################################################################
# Specify location of config-file
##################################################################################
#cd("C://Users//GOVAERTN//OneDrive - VITO//RepresentativeDaysFinder")
cd("C://Users//u0115926//Documents//RepresentativeDaysFinder")


config_file = normpath(joinpath(@__DIR__, "scenarios", "DE_DK_2015_1.yaml"))

findRepresentativeDays(config_file)
