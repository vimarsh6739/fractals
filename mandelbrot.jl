using Plots
exprange(start, stop, len) = exp.(range(log(start), log(stop), length=len))
function mandelbrot(z; lim=75) w = z
    for n = 1:lim;  abs2(w) < 4 ? w = w^2 + z : return n end
    lim + 1
end
x₀, y₀ = -0.5626805, 0.6422555
anim = @animate for (r, l) in zip(exprange(2, 1.35e-6, 120), exprange(100, 2500, 120))
    x = range(x₀-r, x₀+r; length=600); y = range(y₀-r, y₀+r; length=600);
    heatmap(x, y, -log.(log.(mandelbrot.(x' .+ y .* im; lim=round(l))));
        legend=:none, border=:none, ticks=:none, size=(600,600), ratio=1)
end
g = gif(anim; fps=12)
