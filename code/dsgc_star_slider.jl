# SLIDER: Delay master stability function & simulations for the DSGC model (4-node star network)
# Author: Reyk
# last edited: January 27, 2021

# import required packages
using Plots, LaTeXStrings
using Interact, Mux #for slider
using DifferentialEquations, NLsolve
using DelimitedFiles # to read txt files

# slider values (delay in seconds)
tau_list = [0,0.33,0.36,1.7,2.,2.52,2.54,4.5,5,6.85,8.5,8.8,9.9]

# read delay master stability function data
tipping = readdlm("2021/schaefer_tipping.txt", ',', Float64) # boundaries of stability windows
schaefer_sigma = readdlm("2021/schaefer_sigma.txt", ',', Float64) # sigma_max
sch1 = readdlm("2021/schaefer_1.txt", ',', Float64) # eigenvalue lambda_1
sch2 = readdlm("2021/schaefer_2.txt", ',', Float64) # eigenvalue lambda_2
sch4 = readdlm("2021/schaefer_4.txt", ',', Float64) # eigenvalue lambda_4

#model parameters
system_size = 4
P = 1.; K = 8.; alpha = 0.1; gamma = 0.25
power_vector = [3*P, -P, -P, -P]

# adjacency matrix A (star)
adjacency = zeros(system_size,system_size)
for i in 2:system_size
    adjacency[i,1] = K
    adjacency[1,i] = K
end

# define model equation
function four_node_coupling(x, i)
    sum = 0
    for j in 1:system_size
        sum += adjacency[i,j] * sin(x[2*j-1] - x[2*i-1])
    end
    return sum
end

function four_node_coupling2(angles, i)
    sum = 0
    for j in 1:system_size
        sum += adjacency[i,j] * sin(angles[j] - angles[i])
    end
    return sum
end

function four_node_model(dx, x, h, p, t)
    dx[1] = x[2]
    dx[2] = power_vector[1] - p[2]*x[2] + four_node_coupling(x,1) - p[4]*h(p, t-p[end])[2]
    dx[3] = x[4]
    dx[4] = power_vector[2] - p[2]*x[4] + four_node_coupling(x,2) - p[4]*h(p, t-p[end])[4]
    dx[5] = x[6]
    dx[6] = power_vector[3] - p[2]*x[6] + four_node_coupling(x,3) - p[4]*h(p, t-p[end])[6]
    dx[7] = x[8]
    dx[8] = power_vector[4] - p[2]*x[8] + four_node_coupling(x,4) - p[4]*h(p, t-p[end])[8]
end

# fixed point equation
function get_angles!(f,x)
    for i in 1:system_size
        f[i] = four_node_coupling2(x, i) + power_vector[i]
    end
end

# compute phase angles of synchronous fixed point
fix_angles = nlsolve(get_angles!, zeros(4)).zero
fix_angles_shift = zeros(4)
fix_angles_shift .= fix_angles .- fix_angles[1]

# initialize fixed point array
fp = zeros(8)
for i in 1:4
    fp[2*i-1] = fix_angles_shift[i]
end

# random frequency perturbation
for i in 1:system_size
    fp[2i] = 0.2*rand() - 0.1
end

# simulation setup
tspan = (0.0,1000) # number of timesteps (timestep = 0.01 seconds -> duration of simulation = 10 seconds)
x_0 = fp #initial state (random perturbation from fixed point)
history(p,t) = fp # history function (choose random perturbation as history function)

# run simulation for each slider value
mp = @manipulate throttle=.05 for tau in tau_list
    lag2 = [tau]
    p2 = [P, alpha, K, gamma, tau]
    prob2 = DDEProblem(four_node_model, x_0, history, tspan, p2; constant_lags=lag2)
    alg = MethodOfSteps(Tsit5())
    sim2 = solve(prob2, alg, saveat=0.01);
    sim2[1:2:end, :] = mod2pi.(sim2[1:2:end, :] .+ pi) .- pi;

    l = @layout [a b]
    # left plot (delay master stability function)
        plot1 = plot([0,10],[0,0], color="black", linestyle=:dash, lw=1, aspect_ratio=1.6, dpi=200,
            title="Delay master stability function",
            xlabel="delay (seconds)",
            ylabel="dMSF",
            label="",
            yticks=([-1,0,1]),
            ylims=(-2.5,1),
            xlims=(0,10),
            tickfontsize=8,
            titlefontsize=10,
            thickness_scaling=1,
            guidefontsize=9,
            legendfontsize=9,
            legend=:bottomright,
            grid=:none,
            framestyle=:box,
            tickdirection=:out,
            top_margin=0mm,
            right_margin=3mm,
            left_margin=8mm,
            )
            plot!([0,tipping[1]], [0,0], fillrange=(-3,2), fillcolor="black", fillalpha=0.15, linealpha=0, label="")
            plot!([tipping[2],tipping[3]], [0,0], fillrange=(-3,2), fillcolor="black", fillalpha=0.15, linealpha=0, label="")
            plot!([tipping[4],tipping[5]], [0,0], fillrange=(-3,2), fillcolor="black", fillalpha=0.15, linealpha=0, label="")
            plot!([tipping[6],tipping[7]], [0,0], fillrange=(-3,2), fillcolor="black", fillalpha=0.15, linealpha=0, label="")
            plot!(0.001:0.001:10,-sch1, color="navy", ls=:dot, lw=1, label=:none)
            plot!(0.001:0.001:10,-sch2, color="navy", ls=:dot, lw=1, label=:none)
            plot!(0.001:0.001:10,-sch4, color="navy", ls=:dot, lw=1, label=:none)
            plot!(0.001:0.001:10, schaefer_sigma, lw=2, color="navy", label=L"\sigma_{max}(\tau)")
            plot!([tau,tau],[-10,2], lw=2, color="firebrick",label=:none)
    # right plot (simulation)
        plot2 = plot(log2.(0:0.01:1000), sim2[2,:], xlims=(log2(0.09),log2(1000)), ylims=(-0.19,0.19), dpi=200, aspect_ratio=25,
            lw=1,
            label=L"\omega_1",
            color="orange", linealpha=1,
            title="Simulation",
            xlabel="time (seconds)",
            ylabel="frequency deviation (Hz)",
            xticks=(log2.([0.1,1,10,100]), ["0.1","1","10","100"]),
            yticks=([-0.1,0,0.1], ["-0.1","0","0.1"]),
            tickfontsize=8,
            titlefontsize=10,
            thickness_scaling=1,
            guidefontsize=9,
            legendfontsize=8,
            legend=:topright,
            grid=:none,
            framestyle=:box,
            tickdirection=:out,
            top_margin=0mm,
            right_margin=8mm,
            left_margin=3mm,
            )
            plot!(log2.(0:0.01:1000), sim2[4,:], color="dodgerblue", lw=1, label=L"\omega_2")
            plot!(log2.(0:0.01:1000), sim2[6,:], color="firebrick", lw=1, label=L"\omega_3")
            plot!(log2.(0:0.01:1000), sim2[8,:], color="mediumseagreen", lw=1, label=L"\omega_4")

        plot(plot1,plot2, layout=l, size=(800,320))
end

# make slider displayable in webbrowser
ui = dom"div"(mp)
WebIO.webio_serve(page("/", req -> ui), 8001)

# To access the slider, type "localhost:8001/" into your webbrowser
