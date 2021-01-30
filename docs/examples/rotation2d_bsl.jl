# # Rotation of a gaussian distribution
# 
#md # [`notebook`](https://@__NBVIEWER_ROOT_URL__notebooks/rotation_2d_bsl.ipynb)
# 
# ```math
# \frac{df}{dt} +  (y \frac{df}{dx} - x \frac{df}{dy}) = 0
# ```
# 
# 

import SemiLagrangian: advection!, CubicSpline
import VlasovBase: UniformMesh
import SplittingOperators: @Magic
using Plots
pyplot()

#------------------------------------------------------------------------------

function with_bsl(tf::Float64, nt::Int)

   n1, n2 = 32, 64
   mesh1 = UniformMesh(-π, π, n1)
   mesh2 = UniformMesh(-π, π, n2)
   x = mesh1.points
   y = mesh2.points

   dt = tf/nt

   f = zeros(Float64,(n1,n2))

   for (i, xp) in enumerate(x), (j, yp) in enumerate(y)
       xn = cos(tf)*xp - sin(tf)*yp
       yn = sin(tf)*xp + cos(tf)*yp
       f[i,j] = exp(-(xn-1)*(xn-1)/0.2)*exp(-(yn-1)*(yn-1)/0.2)
   end

   anim = @animate for n=1:nt
       
	   @Magic(advection!( f,  mesh1,  y, dt, CubicSpline(), 1),
		    advection!( f,  mesh2, -x, dt, CubicSpline(), 2)
		   )

      surface(f)
                                      
   end

   gif(anim, "rotanim.gif", fps=15); nothing #hide

end

#------------------------------------------------------------------------------

@time f = with_bsl( 2π, 6)


# ![](rotanim.gif)
