abstract type AbstractMethod end

export BSLSpline

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct BSLSpline <: AbstractMethod

    p :: Int

end

"""
$(SIGNATURES)

Return the value at x in [0,1] of the B-spline with integer nodes of degree p with support starting at j.
Implemented recursively using the [De Boor's Algorithm](https://en.wikipedia.org/wiki/De_Boor%27s_algorithm)

```math
B_{i,0}(x) := \\left\\{
\\begin{matrix}
1 & \\mathrm{if}  \\quad t_i ≤ x < t_{i+1} \\\\
0 & \\mathrm{otherwise} 
\\end{matrix}
\\right.
```

```math
B_{i,p}(x) := \\frac{x - t_i}{t_{i+p} - t_i} B_{i,p-1}(x) 
+ \\frac{t_{i+p+1} - x}{t_{i+p+1} - t_{i+1}} B_{i+1,p-1}(x).
```
"""
function bspline(p::Int, j::Int, x::Float64)
   if p == 0
      j == 0 ? (return 1.0) : (return 0.0)
   else
      w = (x - j) / p
      w1 = (x - j - 1) / p
   end
   return w * bspline(p - 1, j, x) + (1 - w1) * bspline(p - 1, j + 1, x)
end


