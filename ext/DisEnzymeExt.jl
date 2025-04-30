module DisEnzyme
import Discontinuity: _gradient
using Enzyme
using Enzyme: Reverse, gradient

Discontinuity._gradient(f, x) = gradient(Reverse, f, x)[1]

end