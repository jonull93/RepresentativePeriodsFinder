############################## Copyright (C) 2019  #############################
#       The content of this file is VITO (Vlaamse Instelling voor              #
#       Technologisch Onderzoek  N.V.) proprietary.                            #
################################################################################
@testset "TimeSeries: DummyColumn" begin
    # Initialize empty time series
    ts = dft.time_series["DummyColumn"]

    @test isa(ts.data, Array{Float64,1})
    @test ts.data == Array{Float64,1}(1:8760)

    @test isa(ts.matrix_full, Array{Float64,2})
    @test ts.matrix_full == reshape(Array{Float64,1}(1:8760), 24, 365)'
#
#     # 79.1667  20.8333   0.0      0.0      0.0
#     # 0.0      58.3333  41.6667   0.0      0.0
#     # 0.0       0.0     37.5     62.5      0.0
#     # 0.0       0.0      0.0     16.6667  83.3333
#     TEST_matrix_bins = [19 5 0 0 0; 0 14 10 0 0; 0 0 9 15 0; 0 0 0 4 20]/24.*100
#     @test isapprox(ts.matrix_bins, TEST_matrix_bins)
#
#     # 79.1667  100.0     100.0  100.0     100.0
#     #  0.0      58.3333  100.0  100.0     100.0
#     #  0.0       0.0      37.5  100.0     100.0
#     #  0.0       0.0       0.0   16.6667  100.0
#     TEST_matrix_bins_cumsum_day = cumsum(TEST_matrix_bins,dims=2)
#     @test isapprox(ts.matrix_bins_cumsum_day, TEST_matrix_bins_cumsum_day)
#
#     WEIGHT = Dict{String, Float64}()
#     RepresentativeDaysFinder.weight!(WEIGHT, ts)
#     @test WEIGHT[ts.name] == 1.
#
#     AREA_TOTAL = Dict{String, Float64}()
#     RepresentativeDaysFinder.area_total!(AREA_TOTAL, ts)
#     @test AREA_TOTAL[ts.name] == 97 * 48.
#
#     AREA_TOTAL_DAYS = Dict{Tuple{String,String}, Float64}()
#     RepresentativeDaysFinder.area_total_days!(AREA_TOTAL_DAYS, dft.periods, ts)
#     @test AREA_TOTAL_DAYS[ts.name, "p001"] == sum(1:24)
#     @test AREA_TOTAL_DAYS[ts.name, "p002"] == sum(25:48)
#     @test AREA_TOTAL_DAYS[ts.name, "p003"] == sum(49:72)
#     @test AREA_TOTAL_DAYS[ts.name, "p004"] == sum(73:96)
#
#     # 19.7917     39.5833     59.375     79.1667     100.0
#     TEST_matrix_bins_cumsum = cumsum(sum([19 5 0 0 0; 0 14 10 0 0; 0 0 9 15 0; 0 0 0 4 20],dims=1)/96.*100,dims=2)[:]
#     @test isapprox(ts.matrix_bins_cumsum, TEST_matrix_bins_cumsum)
#
#     A = Dict{Tuple{String,String,String},Float64}()
#     RepresentativeDaysFinder.cum_bin_end!(A, dft.periods, dft.bins, ts)
#     @test isapprox(A[ts.name, "p001", "b001"],  79.16666666666667)
#     @test isapprox(A[ts.name, "p001", "b002"], 100.)
#     @test isapprox(A[ts.name, "p001", "b003"], 100.)
#     @test isapprox(A[ts.name, "p001", "b004"], 100.)
#     @test isapprox(A[ts.name, "p001", "b005"], 100.)
#
#     @test isapprox(A[ts.name, "p002", "b003"], 100.)
#
#     L = Dict{Tuple{String,String},Float64}()
#     RepresentativeDaysFinder.cum_bin_total!(L, dft.bins, ts)
#     @test isapprox(L[ts.name, "b001"], 19.79166666666667)
#
#     periods = RepresentativeDaysFinder.get_mandatory_periods(ts, dft)
#     @test length(periods) == 2
#     @test periods == String["p001","p004"]
#
# end
#
# ##################################################################################
# # Testing of Diff-methods
# ##################################################################################
# ts_config = config["time_series"][2]
# ts2 = RepresentativeDaysFinder.TimeSeries(dft, ts_config)
#
# ts_diff = RepresentativeDaysFinder.TimeSeries(dft, ts2)
# @testset "TimeSeries: diff1_testCurve_02" begin
#
#     TEST_data = [0.; -47:47]
#     @test isapprox(ts_diff.data, TEST_data)
#
#     periods = RepresentativeDaysFinder.get_mandatory_periods(ts_diff, dft)
#     @test length(periods) == 0
#     @test periods == String[]
#
# end
#
# ##################################################################################
# # Testing of correlation-method
# ##################################################################################
# ts_corr = RepresentativeDaysFinder.TimeSeries(dft, ts, ts2)
#
#
# ##################################################################################
# # Plotting testing
# ##################################################################################
# if false
#     df = DataFrame(x = Float64[],y = Float64[], legend = String[])
#     @show df
#     append!(df, DataFrame(x = range(1,length(ts.data)) / length(ts.data) * 100, y = sort(ts.data, rev = true), legend = "original"))
#
#     dft.u = Dict{String,Int}()
#     dft.u["p001"] = 0
#     dft.u["p002"] = 1
#     dft.u["p003"] = 0
#     dft.u["p004"] = 1
#
#     dft.w = Dict{String,Float64}()
#     dft.w["p001"] = 0.
#     dft.w["p002"] = 2.5
#     dft.w["p003"] = 0.
#     dft.w["p004"] = 1.5
#
#     period_idx = [dft.u[k] for k in sort([k for k in keys(dft.u)])]
#     @show period_idx
#     idx = [i for i in period_idx .* [i for i in 1:length(period_idx)] if i > 0]
#     y = ts.matrix_full[idx, :]'[:]
#     x = ([dft.w[k] for k in sort([k for k in keys(dft.w)]) if dft.w[k] > 0] * ones(1, 24))'[:] / length(ts.data) * 100
#     df2 = sort(DataFrame(x = x, y = y, legend = "reduced"), cols = :y, rev = true)
#     df2[:x] = cumsum(df2[:x])
#     append!(df, df2)
#     p = plot(df, x = :x,y = :y,color = :legend, Geom.line,
#             Guide.xlabel("Duration [-]"), Guide.ylabel("Curve"), Guide.title(ts.name),
#             Coord.Cartesian(xmin = 0,xmax = 100))
#     file_pdf = @sprintf "../results/2017_interconnected_cm_2015/plots/%s.pdf" ts.name
#     draw(PDF(file_pdf, 27cm, 21cm), p)
end
#
#
# # include("testing_real_profiles.jl")
