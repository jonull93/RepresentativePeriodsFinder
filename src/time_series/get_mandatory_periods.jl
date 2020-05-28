############################## Copyright (C) 2019  #############################
#       The content of this file is VITO (Vlaamse Instelling voor              #
#       Technologisch Onderzoek  N.V.) proprietary.                            #
################################################################################
function get_mandatory_periods(ts::TimeSeries, dft::DaysFinderTool)
    periods = []
    if isdefined(ts, :config) && "mandatory_periods" in keys(ts.config)
        for method in ts.config["mandatory_periods"]
            if method == "max"
                val, idx = findmax(ts.matrix_full)
                push!(periods, dft.periods[idx[1]])
            elseif method == "min"
                val, idx = findmin(ts.matrix_full)
                push!(periods, dft.periods[idx[1]])
            end
        end
        @debug("Mandatory periods for $(ts.name): $periods")
    end
    return unique(periods)
end
