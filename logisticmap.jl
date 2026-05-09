using GLMakie, Makie

GLMakie.activate!()

logisticmap(x, R) = R * x * (1 - x)

R  = Observable(2.0)
x0 = Observable(0.2)
x = Observable([x0[]])
cobweb_pts = Observable(Point2f[(x0[], 0)])


fig = Figure(size = (1250, 900))

control_panel = GridLayout()
plot_panel = GridLayout()

fig[1, 1] = control_panel
fig[1, 2] = plot_panel

colsize!(fig.layout, 1, Fixed(200))
colsize!(fig.layout, 2, Auto())

# --- Controls ---

Label(control_panel[1, 1], "Growth Rate (R)", font=:bold)
r_slider = Slider(control_panel[2, 1], range=0:0.01:4, startvalue=2.0)
r_val_label = Label(control_panel[2, 2], @lift string(round($R, digits=2)))

Label(control_panel[3, 1], "Initial Value (x0)", font=:bold)
x0_slider = Slider(control_panel[4, 1], range=0:0.01:1, startvalue=0.2)
x0_val_label = Label(control_panel[4, 2], @lift string(round($x0, digits=2)))

step_btn  = Button(control_panel[5, 1], label="Step")
reset_btn = Button(control_panel[6, 1], label="Reset")
bifdiag_btn = Button(control_panel[7, 1], label="Bifurcation Diagram")

# Sliders, buttons and labels spacing
rowsize!(control_panel, 1, Fixed(20))
rowsize!(control_panel, 2, Fixed(30))
rowsize!(control_panel, 3, Fixed(20))
rowsize!(control_panel, 4, Fixed(30))
rowsize!(control_panel, 5, Fixed(50))
rowsize!(control_panel, 6, Fixed(50))
rowsize!(control_panel, 7, Fixed(50))

rowgap!(control_panel, 10)

# --- Plots ---

ax1 = Axis(plot_panel[1, 1],
    title = "Evolution",
    xlabel = "Iteration",
    ylabel = "x"
)

plotlen = 15
xlims!(ax1, 0, 15)
ylims!(ax1, 0, 1)

ax2 = Axis(plot_panel[2, 1],
    title = "Cobweb Plot",
    xlabel = "x",
    ylabel = "f(x)",
)

xlims!(ax2, 0, 1)
ylims!(ax2, 0, 1)

rowsize!(plot_panel, 1, Fixed(300))
rowsize!(plot_panel, 2, Fixed(300))
rowgap!(plot_panel, 50)

# --- Interaction Logic ---

on(r_slider.value) do val
    R[] = val
end

on(x0_slider.value) do val
    x0[] = val
end

function reset_state!()
    x[] = [x0[]]
    cobweb_pts[] = [Point2f(x0[], 0)] 
    notify(x)
    notify(cobweb_pts)
    xlims!(ax1, 0, 15)
end

onany(R, x0, reset_btn.clicks) do _, _, _
    reset_state!()
end

on(step_btn.clicks) do _
    x_old = last(x[])
    x_new = logisticmap(x_old, R[])

    push!(x[], x_new)

    push!(cobweb_pts[], Point2f(x_old, x_old))
    push!(cobweb_pts[], Point2f(x_old, x_new))

    notify(x)
    notify(cobweb_pts)

    if length(x[]) > plotlen
        xlims!(ax1, 0, length(x[]) + 4)
        global plotlen += 5
    end
end

on(bifdiag_btn.clicks) do _
    R_steps = 10000
    iters_per_R = 250 
    Rs = range(1.0, 4.0, length=R_steps)
    
    plot_Rs = Vector{Float64}(undef, R_steps * iters_per_R)
    plot_xs = Vector{Float64}(undef, R_steps * iters_per_R)
    
    idx = 1
    for (i, r) in enumerate(Rs)
        x = 0.5
        
        for _ in 1:250 # Warmup
            x = r * x * (1 - x) 
        end
        
        for j in 1:iters_per_R # Collect data
            x = r * x * (1 - x)
            plot_Rs[idx] = r
            plot_xs[idx] = x
            idx += 1
        end
    end

    fig = Figure(size=(1000, 600))
    ax = Axis(fig[1, 1], title="Logistic Map Bifurcation", 
                xlabel="R", ylabel="x")
    vals = [1.0, 2.0, 3.0, 3.449, 3.57, 4.0]
    labels = ["1.0", "2.0", "(Period 2) 3.0", "(Period 4) 3.449", "(Chaos) 3.5699", "4.0"]
    ax.xticks = (vals, labels)
    ax.xticklabelrotation = 0.9 
    ax.xticklabelcolor = :black 
    scatter!(ax, plot_Rs, plot_xs, markersize=0.3, color=(:black, 0.25))
    display(fig) 
end

# --- Plotting ---

# Evolution plot
lines!(ax1, x, color=:blue, linewidth=2)

# y = x
xs = 0:0.01:1
lines!(ax2, xs, xs, linestyle=:dash, color=:black)

# x0
scatter!(x0, 0, marker=:circle, color=:blue, markersize=10)

# Cobweb path
lines!(ax2, cobweb_pts, color=:blue, linewidth=1.5)

# Logistic curve 
y_parabola = @lift [logisticmap(val, $R) for val in xs]
lines!(ax2, xs, y_parabola, color=:red, linewidth=2)

display(fig)
wait(display(fig))