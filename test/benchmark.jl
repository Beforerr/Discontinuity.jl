using Discontinuity: norm_std
using Chairmarks

rdata = rand(10000, 3)

@b norm_std(data)


# @b map(nanstd, eachcol(rdata))
f1(v) = norm(nanstd(rdata, dim=1))
f2(v) = norm(nanstd.(eachcol(v)))
f3(v) = norm([nanstd(col) for col in eachcol(v)])
f4(v) = norm(SVector{3}(nanstd(view(v, :, col)) for col in 1:3))

@b f1(rdata), f2(rdata), f3(rdata), f4(rdata)

@b norm(map(nanstd, eachcol(rdata)))

@b detect_variance(data, tau)
@profview_allocs detect_variance(data, tau)
@descend detect_variance(data, tau)