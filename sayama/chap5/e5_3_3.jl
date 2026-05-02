using Plots

gr()

# Model:
# N_t = N_{t-1} + rN_{t-1} (1 - N_{t-1}/K) | r = 1, K = 1, N0 = 0.1
r = 2.5; K = 1
f(x) = x + r*x * (1 - x/K) 

function cobweb(f, x0, steps)
    xs = Float64[]
    ys = Float64[]
    x = x0
    push!(xs, x); push!(ys, x)
    
    for _ in 1:steps
        y = f(x)
        push!(xs, x); push!(ys, y)
        push!(xs, y); push!(ys, y)
        x = y
    end

    return xs, ys
end

x = 0:0.1:1.5

p = plot(x, f.(x), label="Model")
plot!(p, x, x, label="Equilibrium")

x0 = 0.1
steps = 20
xs, ys = cobweb(f, x0, steps)

plot!(p, xs, ys, label="Cobweb", linestyle=:dash)

display(p)

println("Press Enter...")
readline()