using Documenter
using RepresentativePeriodsFinder

makedocs(
    sitename = "RepresentativePeriodsFinder",
    format = Documenter.HTML(),
    modules = [RepresentativePeriodsFinder],
    pages = [
        "Home" => "index.md",
        "User Guide" => [
            "Getting started" => "getting_started.md",
            "Loading time series" => "loading_time_series.md",
            "Selecting periods" => "selecting_periods.md",
            "Outputs" => "outputs.md",
        ],
        "Methods" => "methods.md",
        "API Reference" => "api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
