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
    # df = DataFrame(OrderedDict(:Timestep => Int32[], colName => Union{Missing,Float64}[]))
    df = DataFrame([Int,Union{Missing,Float64}], [:Timestep,colName])
    nt = size(data)[1]
    for i = 1:nt
        try
            if typeof(data[i,colNum]) == String
                push!(df, [i parse(Float64, data[i,colNum])])
            elseif typeof(data[i,colNum]) <: Number
                push!(df, [i Float64(data[i,colNum])])
            else
                println("Probably wrong")
            end
        catch
            push!(df, [i missing])
        end
    end
    return df
end

"""
    Fills missing rows on temporal axis of dataframe in datecolumn.
    resultion in hours
"""
function check_temporal_consistency(df::DataFrame, resultion::Float64, format::String, datecolumnOrig:: Symbol; datecolumn = :DateFormated)
    dformat = DateFormat(format)
    # df_str = "y-m-d HH:MM:SSz"

    # put all date in same time zone and add column to df
    dt_df = Dates.Minute(resultion*60.0)
    tzone = tz"Europe/Brussels"
    df[datecolumn] = map(str->astimezone(ZonedDateTime(str, dformat), tzone),df[datecolumnOrig])

    # create empty row df
    @show empty_row = similar(df, 1)
    for n in names(empty_row)
        if n != datecolumn
            empty_row[n] = missing
        end
    end

    #run through df and add rows
    for i in 2:length(df[datecolumn])
        t_current = df[datecolumn][i]
        t_prev = df[datecolumn][i-1]
        dt =  Dates.Minute(t_current -  t_prev)
        if  dt != dt_df
            println((i, df[datecolumn][i], dt))
            #add line
            closed = false
            t_add = t_prev
            while !closed
                tba = deepcopy(empty_row)
                t_add = t_add + dt_df
                tba[datecolumn] = t_add
                append!(df, tba)
                if Dates.Minute(t_current - t_add) == dt_df
                    closed = true
                end
            end
        end
    end

    sorted_df = sort(df, datecolumn)
    return sorted_df
end

"""
    Fills in any missing values in the data frame using linear interpolation.
"""
function interpolatedataframe(df::DataFrame, colName::Symbol)

    itp = makeinterpolator(df)

    # Make the data frame without missing values
    @show nt = size(df)[1]
    df = DataFrame([[Int32(i) for i in 1:nt], [itp(i) for i in 1:nt]], [:Timestep, colName])

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
    df = DataFrame([[Int32(i) for i in 1:nt*ts1/ts2], [itp(i) for i in 1:ts2/ts1:nt+1-ts2/ts1]], [:Timestep, colName])

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
