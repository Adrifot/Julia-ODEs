using Plots

gr()

# Model 1:
# x_t = x_{t-1} + 0.1 | x0 = 0.1
f1(x) = x + 0.1

# Model 2:
# x_t = 1.1x_{t-1} | x0 = 0.1
f2(x) = 1.1x

function cobweb(f, x0, n)
    xs = Float64[]
    ys = Float64[]
    x = x0
    push!(xs, x); push!(ys, x)
    
    for _ in 1:n
        y = f(x)
        push!(xs, x); push!(ys, y)
        push!(xs, y); push!(ys, y)
        x = y
    end

    return xs, ys
end

x = 0:0.1:2

p = plot(x, f1.(x), label="Model 1")
p = plot!(x, f2.(x), label="Model 2")
plot!(p, x, x, label="Equilibrium")

x0 = 0.1
steps = 20
xs1, ys1 = cobweb(f1, x0, steps)
xs2, ys2 = cobweb(f2, x0, steps+10)

plot!(p, xs1, ys1, label="M1 Cobweb", linestyle=:dash)
plot!(p, xs2, ys2, label="M2 Cobweb", linestyle=:dash)

display(p)

println("Press Enter...")
readline()