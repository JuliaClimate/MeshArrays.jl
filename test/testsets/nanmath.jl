@testset "nanmath" begin
    x=[NaN 1 2]
    nansum(x)
    nansum(x,1)
    nanmax(x,2)
    nanmin(x,2)

    nanmean(NaN,1)
    nanmean(1,NaN)
    nanmean(2,1)
    nanmean(NaN,NaN)
end
