# Fractal animations with Enzyme

This repository contains two complex-plane zoom animations:

- Mandelbrot adapted from a [Julia community example](https://discourse.julialang.org/t/seven-lines-of-julia-examples-sought/50416/47)
- Newton-fractal implementation that uses Enzyme to differentiate
  the selected polynomial.

Both animations evaluate an orbit for every pixel, color the result by the
number of iterations required to escape or converge, and progressively shrink
the complex-plane viewport.

## Setup

Instantiate the Julia environment from the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Plots` renders the heatmaps and GIFs. The Newton implementation additionally
uses `Enzyme` for automatic differentiation.

## Mandelbrot animation

[`mandelbrot.jl`](mandelbrot.jl) samples each pixel as the complex parameter
`c` in

```math
z_{n+1} = z_n^2 + c.
```

The escape-time function returns the first iteration for which
`abs2(z) >= 4`. The animation zooms toward
`-0.5626805 + 0.6422555im`, increasing the iteration limit as the viewport gets
smaller so that the boundary remains detailed.

Run it with:

```sh
julia --project=. mandelbrot.jl
```

## Newton-fractal animation

[`newton.jl`](newton.jl) visualizes Newton iteration for the polynomial

```math
p(z) = z^8 + 15z^4 - 16.
```

For every pixel, the pixel coordinate supplies an independent starting value
`z_0`. The orbit follows

```math
z_{n+1} = z_n - \frac{p(z_n)}{p'(z_n)}
```

until `abs(p(z))` falls below the convergence tolerance or the iteration limit
is reached. Enzyme computes the complex derivative with forward-mode automatic
differentiation. Because `p` is holomorphic, seeding the derivative with
`1 + 0im` produces the usual complex derivative `p'(z)`.

The heatmap stores the convergence iteration for every pixel. `log1p` is
applied only to spread those integer counts more evenly across the color
palette.

Run the animation with Julia threads enabled:

```sh
julia --threads=auto --project=. newton.jl
open newton_fractal.gif
```

The default animation renders 120 frames at 600 by 600 pixels and zooms around
a repelling period-2 point of the Newton map:

```julia
0.5309620407446378 + 0.19807575472155833im
```

Unlike an attracting root, a repelling periodic point lies on the Julia set,
so arbitrarily small neighborhoods continue to intersect different attracting
basins.

This is what the final image looks like:

![Newton Fractal](./newton_fractal.gif)

## Finding 'good' Newton zoom centers

Ironically, finding a good zoom center was harder than I thought. I initially spammed 
center values, then saw how the animation turned out. A little bit of digging had me 
solve the equation $N(z) = 0$, where $N(z)$ is the Newton-Raphson iteration expressed. 

[`find_period2_centers.jl`](find_period2_centers.jl) reproduces the search used
to choose the default Newton zoom center. It numerically solves

```math
N^2(z) - z = 0,
```

where `N` is the Newton map. It then removes the attracting fixed points,
deduplicates numerical solutions, and retains cycles whose multiplier satisfies

```math
\left|(N^2)'(z)\right| > 1.
```

Print the selected center and its cycle information with:

```sh
julia --project=. find_period2_centers.jl
```

Print every candidate found by the default grid search with:

```sh
julia --project=. find_period2_centers.jl --all
```

The final selection heuristic favors a strongly repelling point away from the
real and diagonal symmetry axes. Other reported candidates can be passed as
the `center` keyword to `render_newton_animation`.
