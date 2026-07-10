"""
Search for repelling period-2 points of the Newton map associated with

    p(z) = z^8 + 15z^4 - 16.

Repelling periodic points lie on the Julia set, so they are exact candidates
for a persistent zoom center. Run with `--all` to print every point found by
the grid search.
"""

p(z) = z^8 + 15z^4 - 16
dp(z) = 8z^7 + 60z^3
ddp(z) = 56z^6 + 180z^2

newton_map(z) = z - p(z) / dp(z)

# Differentiating N(z) = z - p(z)/p'(z) gives this compact expression.
newton_map_derivative(z) = p(z) * ddp(z) / dp(z)^2

struct PeriodicPoint
    point::ComplexF64
    multiplier::ComplexF64
end

"""Return `N^period(z)` and the derivative `(N^period)'(z)."""
function iterate_with_derivative(z, period::Integer)
    value = z
    multiplier = one(z)

    for _ in 1:period
        multiplier *= newton_map_derivative(value)
        value = newton_map(value)

        if !isfinite(value) || !isfinite(multiplier)
            return nothing
        end
    end

    return (; value, multiplier)
end

"""
Solve `N^period(z) - z = 0` from one initial seed with Newton's method.
The derivative of that equation is `(N^period)'(z) - 1`.
"""
function solve_periodic_point(
    seed,
    period::Integer;
    tolerance::Real=1e-12,
    derivative_tolerance::Real=1e-14,
    limit::Integer=100,
)
    z = ComplexF64(seed)

    for _ in 1:limit
        result = iterate_with_derivative(z, period)
        result === nothing && return nothing

        residual = result.value - z
        abs(residual) <= tolerance &&
            return PeriodicPoint(z, result.multiplier)

        derivative = result.multiplier - 1
        abs(derivative) > derivative_tolerance || return nothing

        z -= residual / derivative
        isfinite(z) || return nothing
    end

    return nothing
end

"""
Search a square grid of seeds for repelling points of exact period two.

The fixed-point check removes the roots of `p`, while the multiplier check
keeps only repelling cycles. Numerically repeated solutions are deduplicated.
"""
function find_period2_points(;
    lower::Real=-2.5,
    upper::Real=2.5,
    samples::Integer=61,
    fixed_tolerance::Real=1e-7,
    duplicate_tolerance::Real=1e-7,
)
    points = PeriodicPoint[]
    grid = range(lower, upper; length=samples)

    for x in grid, y in grid
        candidate = solve_periodic_point(complex(x, y), 2)
        candidate === nothing && continue

        z = candidate.point

        # N(z) == z would be a period-1 point: one of the attracting roots.
        abs(newton_map(z) - z) > fixed_tolerance || continue

        # A repelling cycle has |(N^2)'(z)| > 1 and lies on the Julia set.
        abs(candidate.multiplier) > 1 + fixed_tolerance || continue

        all(abs(z - existing.point) > duplicate_tolerance for existing in points) ||
            continue

        push!(points, candidate)
    end

    sort!(points; by=candidate -> (angle(candidate.point), abs(candidate.point)))
    return points
end

"""
Reproduce the center selected for the animation preview.

The aesthetic heuristic chooses a point in the first octant, away from the
real and diagonal symmetry axes, with a strong multiplier, then takes the one
nearest the origin. Other candidates may produce equally valid zooms.
"""
function choose_preview_center(points)
    eligible = filter(points) do candidate
        z = candidate.point
        0 < imag(z) < real(z) && abs(candidate.multiplier) > 10
    end

    isempty(eligible) && error("The search did not find a preview candidate")
    sort!(eligible; by=candidate -> abs(candidate.point))
    return first(eligible)
end

function main()
    points = find_period2_points()
    selected = choose_preview_center(points)
    partner = newton_map(selected.point)
    cycle_error = abs(newton_map(partner) - selected.point)

    println("Found $(length(points)) repelling period-2 points")
    println("Selected center:  ", selected.point)
    println("Cycle partner:    ", partner)
    println("Cycle error:      ", cycle_error)
    println("|multiplier|:     ", abs(selected.multiplier))

    if "--all" in ARGS
        println("\nAll candidates:")
        for candidate in points
            println(candidate.point, "    |multiplier| = ", abs(candidate.multiplier))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
