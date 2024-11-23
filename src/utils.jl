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

"""Return the angle between two vectors."""
function vector_angle(v1, v2)
    cos_angle = v1 ⋅ v2 / (norm(v1) * norm(v2))
    return acosd(cos_angle)
end