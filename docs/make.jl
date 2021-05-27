using Documenter, Literate, RepresentativePeriodsFinder

# This is heavily copied from COSMO.jl: https://github.com/oxfordcontrol/COSMO.jl/blob/master/docs/make.jl#L1

@info "Building examples..."
fix_math_md(content) = replace(content, r"\$\$(.*?)\$\$"s => s"```math\1```")
fix_suffix(filename) = replace(filename, ".jl" => ".md")
function postprocess(cont)
      """
      The source files for all examples can be found in [/examples](https://gitlab.kuleuven.be/UCM/representativeperiodsfinder.jl/-/tree/master/examples/).
      """ * cont
end

# Copy some of the tests and input data to the example folder
data_path = joinpath(@__DIR__, "..", "test", "input_data")
test_path = joinpath(@__DIR__, "..", "test")
example_file_names = ["test_days_re_ordering.jl"]
examples_path = joinpath(@__DIR__, "..", "examples")
isdir(examples_path) == false && mkdir(examples_path)
cp(data_path, joinpath(examples_path, "input_data"); force=true)
for ex in example_file_names
    cp(joinpath(test_path, ex), joinpath(examples_path, ex); force=true)
end

# Build the example documentation
build_path =  joinpath(@__DIR__, "src", "examples/")
for file in readdir(examples_path)
    endswith(file, ".jl") == false && continue
    Literate.markdown(
        joinpath(examples_path, file), build_path; preprocess=fix_math_md, 
        postprocess=postprocess, documenter=true, credit=true
    )
end

@info "Building documentation..."
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
