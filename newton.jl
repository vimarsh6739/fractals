using Enzyme
using Plots

# The polynomial whose roots define the attracting basins.
p(z) = z^8 + 15*z^4 - 16*one(z)

"""
    polynomial_derivative(z)

Differentiate `p` at the complex point `z` with Enzyme.  Since `p` is
holomorphic, its directional derivative along `1 + 0im` is the usual complex
derivative p'(z).
"""
@inline function polynomial_derivative(z::Complex{T}) where {T<:AbstractFloat}
    Enzyme.autodiff(
        Enzyme.Forward,
        p,
        Enzyme.Duplicated(z, one(z)),
    )[1]
end

"""
    newton_iterations(z0; limit=75, tolerance=1e-10)

Return the number of Newton steps needed for the orbit beginning at `z0` to
reach a root of `p`.  `limit + 1` means that the orbit did not converge within
the iteration limit.
"""
function newton_iterations(
    z0::Complex{T};
    limit::Integer=75,
    tolerance::Real=1e-10,
) where {T<:AbstractFloat}
    z = z0
    tolerance_squared = T(tolerance)^2

    for iteration in 0:limit
        pz = p(z)
        isfinite(pz) || return limit + 1
        abs2(pz) <= tolerance_squared && return iteration
        iteration == limit && break

        derivative = polynomial_derivative(z)
        (!isfinite(derivative) || iszero(derivative)) && return limit + 1

        z -= pz / derivative
        isfinite(z) || return limit + 1
    end

    return limit + 1
end

"""Evaluate the Newton convergence time for every complex point in a viewport."""
function convergence_image(xs, ys; limit::Integer=75, tolerance::Real=1e-10)
    iterations = Matrix{Int}(undef, length(ys), length(xs))

    # Each pixel supplies its own starting point z0 = x + iy.
    Threads.@threads for row in eachindex(ys)
        y = ys[row]
        @inbounds for column in eachindex(xs)
            z0 = complex(xs[column], y)
            iterations[row, column] = newton_iterations(
                z0;
                limit=limit,
                tolerance=tolerance,
            )
        end
    end

    return iterations
end

geometric_range(start, stop, length) =
    exp.(range(log(start), log(stop); length=length))

a=((45 + sqrt(1577.0)) / 14)^(1/4) / sqrt(2)
"""
    render_newton_animation(; kwargs...)

Render a zoom into the pole at `z = 0`.  The centre of the viewport is not a
single Newton starting point: every pixel is a different starting point.
"""
function render_newton_animation(;
    center::Complex=0.5309620407446378 + 0.19807575472155833im,
    start_radius::Real=0.5,
    end_radius::Real=2e-8,
    frames::Integer=120,
    resolution::Integer=600,
    first_limit::Integer=80,
    final_limit::Integer=250,
    tolerance::Real=1e-10,
    fps::Integer=12,
    filename::AbstractString="newton_fractal.gif",
)
    radii = geometric_range(start_radius, end_radius, frames)
    limits = round.(Int, geometric_range(first_limit, final_limit, frames))

    # Compile the Enzyme derivative before entering a potentially threaded loop.
    polynomial_derivative(0.5 + 0.5im)

    animation = @animate for (radius, limit) in zip(radii, limits)
        xs = range(real(center) - radius, real(center) + radius; length=resolution)
        ys = range(imag(center) - radius, imag(center) + radius; length=resolution)
        iterations = convergence_image(xs, ys; limit=limit, tolerance=tolerance)

        # log1p only spaces the colours more evenly; each value still represents
        # the integer number of iterations required for convergence.
        heatmap(xs,ys,log1p.(iterations);
            color=:inferno,
            legend=:none,
            colorbar=false,
            border=:none,
            ticks=:none,
            size=(resolution, resolution),
            ratio=1,
        )
    end

    return gif(animation, filename; fps=fps)
end

if abspath(PROGRAM_FILE) == @__FILE__
    render_newton_animation()
end
