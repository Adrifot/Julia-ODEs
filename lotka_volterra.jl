using Plots
using Plots.Measures
using DifferentialEquations
gr()

params = 1, 0.04, 0.01, 0.8
dt = 0.01
simtime = 100.0
p0 = 40.0, 3.0

function lotka_volterra_euler(params, p0, dt, simtime)
    α, β, γ, δ = params
    nsteps = Int(simtime / dt)
    
    # Pre-allocation
    ox = collect(0:dt:simtime)
    oy1 = zeros(length(ox))
    oy2 = zeros(length(ox))

    oy1[1], oy2[1] = p0

    for i in 1:nsteps
        x, y = oy1[i], oy2[i]
        
        dx = (α*x - β*x*y) * dt 
        dy = (γ*x*y - δ*y) * dt
        
        next_x = x + dx
        next_y = y + dy
        
        oy1[i+1] = max(0.0, next_x)
        oy2[i+1] = max(0.0, next_y)
    end

    return ox, oy1, oy2
end

function lotka_volterra(du, u, p, t)
    x, y = u 
    α, β, γ, δ = p 
    du[1] = α*x - β*x*y
    du[2] = γ*x*y - δ*y
end

ox, oy1, oy2 = lotka_volterra_euler(params, p0, dt, simtime)

p1 = plot(
    ox, [oy1 oy2], 
    label=["Prey (x)" "Predator (y)"], 
    title="Population vs Time",
    xlabel="Time", ylabel="Population"
)

p2 = plot(
    oy1, oy2, 
    title = "Phase Portrait",
    xlabel = "x (Prey)",
    ylabel = "y (Predator)",
    label = "Orbit"
)
          
scatter!(p2, [oy1[1]], [oy2[1]], label="Start")

combined_plot = plot(
    p1, p2, 
    layout = (1, 2),
    size=(1200, 500),
    margin = 5mm,      
    left_margin = 10mm
)

display(combined_plot)

println("Press Enter...")
readline()

prob = ODEProblem(lotka_volterra, collect(p0), simtime, collect(params))
sol = solve(prob)

p3 = plot(sol, title="Lotka-Volterra Solution", xlabel="Time", ylabel="Population",
            label=["Prey (x)" "Predator (y)"])

p4 = plot(sol, vars=(1, 2), title="Phase Portrait", xlabel="Prey", ylabel="Predator",
            label = "Orbit")

combined_plot2 = plot(
    p3, p4, 
    layout = (1, 2),
    size=(1200, 500),
    margin = 5mm,      
    left_margin = 10mm
)

display(combined_plot2)

println("Press Enter...")
readline()