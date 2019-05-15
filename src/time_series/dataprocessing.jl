# using CSV
# using DataFrames
# using Interpolations


"""
    Reads a .csv file downloaded from the ENTSO-E transparency website into a data frame.
"""
function ENTSOEcsv2dataframe(csvName::String, colNum::Int, colName::Symbol; delim = ',')
    # Read data
    data = CSV.read(csvName, delim = delim)
    # Values are processed as strings, hence the next step
    df = DataFrame(Dict(:Timestep => Int32[], colName => Union{Missing,Float64}[]))
    # df = DataFrame([Int,Union{Missing,Float64}], [:Timestep,colName])
    nt = size(data)[1]
    for i = 1:nt
        try
            push!(df, [i Float64(data[i,colNum])])
        catch
            push!(df, [i missing])
        end
    end
    return df
end

function makeinterpolator(df::DataFrame)
    # make array and interpolation grid filtering out missing values
    y = [j for j in df[:,2] if !ismissing(j)]
    grid = [i for (i,j) in zip(df[:,1], df[:,2]) if !ismissing(j)]

    # Make interpolation object
    return itp = LinearInterpolation(grid, y, extrapolation_bc=Flat())
end

"""
    Fills in any missing values in the data frame using linear interpolation.
"""
function interpolatedataframe(df::DataFrame, colName::Symbol)

    itp = makeinterpolator(df)

    # Make the data frame without missing values
    @show nt = size(df)[1]
    df = DataFrame([[i for i in 1:nt], [itp(i) for i in 1:nt]], [:Timestep, colName])

    # Determine if anything missing, if not return the data frame
    if all(describe(df)[:,:nmissing] .== nothing)
        return df
    else
        error("Interpolation of missing values failed")
    end

end

"""
    Fills in any missing values using linear interpolation and also converts data with a sampling time of ts1 to ts2.
"""
function interpolatedataframe(df::DataFrame, colName::Symbol, ts1, ts2)
    itp = makeinterpolator(df)

    # Make the data frame with interpolated missing values and a new sampling time
    nt = size(df)[1]
    df = DataFrame([[i for i in 1:nt*ts1/ts2], [itp(i) for i in 1:ts2/ts1:nt+1-ts2/ts1]], [:Timestep, colName])

    # Determine if anything missing, if not return the data frame
    if all(describe(df)[:,:nmissing] .== nothing)
        return df
    else
        error("Interpolation of missing values failed")
    end
end

"""
    normalisegeneration(df::DataFrame, colName::Symbol, tg::Number)

    Normalises RES generation profile `df[:,:colName]` by the total available generation capacity `tg`. If any of the resulting values 1, then the profile is again scaled so that all values are <= 1
"""
function normalisegeneration(df::DataFrame, colName::Symbol, tg::Number)
    df[:,colName] = df[:,colName]/tg
    if maximum(df[:,colName]) > 1
        df[:,colName] = 1/maximum(df[:,colName])*df[:,colName]
    end
end
