using Plots
gr()

# x_t = ax_{t-1} + b
# x_0 = 1

function simeq(x0, a, b, tmax)
    x = Vector{Float64}(undef, tmax)
    x[1] = x0
    @inbounds for t in 2:tmax
        x[t] = a * x[t-1] + b
    end
    return x
end

function run(x0, param_pairs, tmax)
    res = Vector{Vector{Float64}}(undef, length(param_pairs))
    l = length(param_pairs)
    for pair_idx in 1:l
        res[pair_idx] = simeq(x0, first(param_pairs[pair_idx]), last(param_pairs[pair_idx]), tmax)
    end
    return res
end

parampairs = [
    1.1 => 0.5,
    0.8 => -0.8,
    1.0 => 1.0,
    -0.5 => 0.25,
    -1.0 => 0.0
]

results = run(1.0, parampairs, 20)

p = plot(title="Dynamics of x_t = ax_{t-1} + b", xlabel="Time (t)", ylabel="x", legend=:topleft)

for (i, pair) in enumerate(parampairs)
    label = "a=$(first(pair)), b=$(last(pair))"
    plot!(p, results[i], label=label, marker=:circle)
end

display(p)

println("Press Enter...")
readline()