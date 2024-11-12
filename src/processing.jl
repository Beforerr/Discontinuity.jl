
"""
Scalar projection, also known as the scalar resolute.

https://en.wikipedia.org/wiki/Scalar_projection
https://en.wikipedia.org/wiki/Vector_projection
"""
sproj(a, b) = dot(a, b) / norm(b)
