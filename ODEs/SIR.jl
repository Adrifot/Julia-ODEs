using Plots
using Plots.Measures
using DifferentialEquations
gr()

"""
SIR Infection model

Equations:
dS/dt = -βSI / N 
dI/dt = βSI / N - γI
dR/dt = γI

Params:
β - infection rate
γ - recovery rate
"""

params = 0.55, 0.15
dt = 0.01
simtime = 50.0
p0 = 999.0, 1.0, 0.0
N = sum(p0)

function sir_euler(params, p0, dt, simtime)
    β, γ = params
    nsteps = Int(simtime/dt)

    ox = collect(0:dt:simtime)
    os = zeros(length(ox))
    oi = zeros(length(ox))
    or = zeros(length(ox))

    os[1], oi[1], or[1] = p0 
    

    for t in 1:nsteps
        S, I, R = os[t], oi[t], or[t]

        dS = (-β*S*I / N) * dt
        dI = (β*S*I/N - γ*I) * dt
        dR = (γ*I) * dt

        nS = S + dS
        nI = I + dI
        nR = R + dR

        os[t+1] = min(max(0.0, nS), N)
        oi[t+1] = min(max(0.0, nI), N)
        or[t+1] = min(max(0.0, nR), N)
    end
    return ox, os, oi, or
end

function sir(du, u, p, t)
    S, I, R = u
    β, γ, N = p
    du[1] = -β*S*I / N
    du[2] = β*S*I / N - γ*I
    du[3] = γ*I 
end

ox, os, oi, or = sir_euler(params, p0, dt, simtime)
sol = ODEProblem(sir, collect(p0), simtime, [params..., N]) |> solve

p1 = plot(
    ox, [os oi or], 
    label=["S" "I" "R"], 
    title="Population vs Time",
    xlabel="Time", ylabel="Population"
)

p2 = plot(
    os, oi, or,
    title = "Phase Portrait",
    xlabel = "S",
    ylabel = "I",
    zlabel = "R",
    label = "Orbit"
)
          
scatter!(p2, [os[1]], [oi[1]], [or[1]], label="Start")

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

p3 = plot(sol, title="SIR Solution", xlabel="Time", ylabel="Population",
            label=["S" "I" "R"])

p4 = plot(sol, vars=(1, 2, 3), title="Phase Portrait", xlabel="S", ylabel="I", zlabel="R",
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