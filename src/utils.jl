const RSQUARED = "fit.stat.rsquared"

"""
    keep_good_fit!(df; rsquared=0.9)

Keep only the rows with a fit statistic R² greater than the given value.
"""
keep_good_fit!(df; rsquared=0.9) = filter!(RSQUARED => >(rsquared), df)
keep_good_fit(df; rsquared=0.9) = filter(RSQUARED => >(rsquared), df)


"""
Scalar projection, also known as the scalar resolute.

https://en.wikipedia.org/wiki/Scalar_projection
https://en.wikipedia.org/wiki/Vector_projection
"""
sproj(v1, v2) = v1 ⋅ v2 / norm(v2)

"""The angle between two vectors"""
function angle_between(v1, v2)
    return acosd(v1 ⋅ v2 / (norm(v1) * norm(v2)))
end

function angle_between_90(v1, v2)
    θ = angle_between(v1, v2)
    ifelse(θ > 90, 180 - θ, θ)
end

const vector_angle = angle_between