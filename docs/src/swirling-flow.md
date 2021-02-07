# Swirling flow

In this problem we use a time-dependent velocity field

```math
\begin{aligned}
  u_x(x,y,t) &= \sin^2(\pi x) \sin(2 \pi y) g(t) \\
  u_y(x,y,t) &= -\sin^2(\pi y) \sin(2 \pi x) g(t)
\end{aligned}
```

This represents a swirling flow that distorts the vorticity field,
reaching a maximum distortion at :math:`t=T/2`. At that point the flow
reverses and the vorticity profile returns to its initial value.

[ref: Gkeyll Simulation Journal](http://ammar-hakim.org/sj/je/je16/je16-ldg.html)
