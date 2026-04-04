using Plots
using Plots.Measures
gr()

α, β, γ, δ = 1, 0.04, 0.01, 0.8
dt = 0.01
simtime = 100.0
nsteps = Int(simtime / dt)

x, y = 40.0, 3.0
pophistory = [(x, y)] 

for t in 1:nsteps

    dx = (α*x - β*x*y) * dt 
    dy = (γ*x*y - δ*y) * dt
    
    global x += dx
    global y += dy
    
    push!(pophistory, (x < 0 ? 0 : x, y < 0 ? 0 : y))
end


ox = range(0, step=dt, length=length(pophistory))
oy1 = first.(pophistory)
oy2 = last.(pophistory)

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
