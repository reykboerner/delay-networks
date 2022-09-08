"""
Script to calculate the delay master stability function (dMSF)
for the Decentral Smart Grid Control (DSGC) model on a star topology

More info: https://github.com/reykboerner/delay-networks
Author: Reyk Börner
Contact: reyk.boerner@reading.ac.uk
"""

# Import packages
using Roots, NLsolve, LinearAlgebra
using Plots, LaTeXStrings, DelimitedFiles

# Setup system
system_size = 4                         # number of nodes
power_vector = [3., -1., -1., -1.]      # net power production/consumption at nodes
start_fp = zeros(4)                     # initial guess for steady state

# Adjacency matrix for star topology
function adj(K)
    A = zeros(system_size,system_size)
    for i in 2:system_size
        A[i,1] = K
        A[1,i] = K
    end
    return A
end

# Coupling term
function coupling_term(angles, i, K)
    sum = 0
    for j in 1:system_size
        sum = sum + adj(K)[i,j] * sin(angles[j] - angles[i])
    end
    return sum
end

# Calculate eigenvalues of Laplacian matrix
function get_laplacian_eigvals(K)
    L = zeros(system_size, system_size)
    strength = K
    function get_angles!(f,x)
        for i in 1:system_size
            f[i] = coupling_term(x,i,strength) + power_vector[i]
        end
    end

    fix_angles = nlsolve(get_angles!, start_fp).zero

    for i in 1:system_size
        rowsum = 0.
        for j in 1:system_size
            if i != j
                L[i,j] = - adj(K)[i,j] * cos(fix_angles[j] - fix_angles[i])
                rowsum = rowsum + L[i,j]
            end
        L[i,i] = - rowsum
        end
    end

    eigenvalues = eigvals(L)
    return eigenvalues
end

# Find decisive root a* for each k = 1...4
function find_a_odd(a,g,m,tau,K,lambdak)
    """
    Finds the decisive root of the characteristic equation for given tau.
        a, g, m: Model parameters α, γ, and inertia m
        tau: delay τ
        K: coupling strength
        lambdak: value of k-th Laplacian eigenvalue λ_k
    """
    q = lambdak * tau^2
    p = a * tau
    # real part of characteristic function
    F(y) = (-y^2 + q) * cos(y) - p*y * sin(y)

    # find all zeros in interval of interest
    a_odd_upperbound = (sqrt(K) + pi) * 10
    F_roots = find_zeros(F, 0., a_odd_upperbound)

    # choose the odd root closest to sqrt_q
    F_roots_odd_dist = zeros(div(length(F_roots), 2))
    for i in 1:length(F_roots_odd_dist)
        F_roots_odd_dist[i] = abs(sqrt(q) - F_roots[2*i])
    end
    a_odd = F_roots[2*argmin(F_roots_odd_dist)]

    return a_odd
end

# Compute stability conditions for each k, for all tau in range
function stability_criterion(a,g,m,K,k,tau_range=0.001:0.001:10)
    """
    Calculates the dMSF for given k in the specified range of delay tau.
        a, g, m: Model parameters α, γ, and inertia m
        K: coupling strength
        k: index of k-th Laplacian eigenvalue λ_k
        tau_range: array of τ values for which to compute the dMSF
    """
    strength = K
    laplacian_eigval = get_laplacian_eigvals(strength)[k]
    a_odd_list = []
    stabil_list = []
    for tau in 0.001:0.001:10
        a_odd_tau = find_a_odd(a,g,m,tau,K,laplacian_eigval)
        stabil_tau = ((a_odd_tau^2 - laplacian_eigval*tau^2)^2 + a^2 * a_odd_tau^2 * tau^2)^(1/2) / a_odd_tau /g/tau - 1
        push!(a_odd_list, a_odd_tau)
        push!(stabil_list, stabil_tau)
    end
    return stabil_list
end

# For each k, execute the function stability_criterion to compute the dMSF
# (Note that sch2 = sch3 due to the symmetry of the star topology)
sch1 = stability_criterion(0.1,0.25,1.,8.,1)
sch2 = stability_criterion(0.1,0.25,1.,8.,2)
sch3 = stability_criterion(0.1,0.25,1.,8.,3)
sch4 = stability_criterion(0.1,0.25,1.,8.,4)

# Find critical tau values where the dmSF changes sign
function get_tc(array)
    intersec = []
    for i in 2:length(array)
        if array[i-1]*array[i] <= 0
            push!(intersec, i-1)
        end
    end
    return float(intersec * 0.001)
end

tc_0 = get_tc(sch1)
tc_1 = get_tc(sch3)
tc_N = get_tc(sch4)
tipping = [tc_N[1], tc_N[4], tc_N[5], tc_N[8], tc_N[9], tc_N[12], tc_N[13], tc_N[16], tc_N[17]]

# Compute supremum of the dMSFs sch1, sch2, sch4
schaefer_sigma = zeros(10000)
for i in 1:10000
    temp = [-sch1[i],-sch2[i],-sch4[i]]
    schaefer_sigma[i] = maximum(temp)
end

# Plot results
plot([0,10],[0,0], color="black", linestyle=:dot, lw=2, aspect_ratio=1.3, dpi=200,
    #xlabel=L"\tau / s",
    #ylabel=L"\sigma(\tau)",
    label="",
    yticks=([-1,0,1]),
    ylims=(-2.5,1),
    xlims=(0,10),
    #ytickfontsize=12,
    tickfontsize=13,
    thickness_scaling=1,
    guidefontsize=16,
    legendfontsize=10,
    legend=:bottomright,
    grid=:true,
    framestyle=:box,
    tickdirection=:in,
    )
    plot!([0,tipping[1]], [0,0], fillrange=(-3,2), fillcolor="black", fillalpha=0.15, linealpha=0, label="")
    plot!([tipping[2],tipping[3]], [0,0], fillrange=(-3,2), fillcolor="black", fillalpha=0.15, linealpha=0, label="")
    plot!([tipping[4],tipping[5]], [0,0], fillrange=(-3,2), fillcolor="black", fillalpha=0.15, linealpha=0, label="")
    plot!([tipping[6],tipping[7]], [0,0], fillrange=(-3,2), fillcolor="black", fillalpha=0.15, linealpha=0, label="")
    plot!(0.001:0.001:10,-sch1, color="orange", lw=2, label=L"\sigma(\lambda_1, \tau)")
    plot!(0.001:0.001:10,-sch3, color="mediumseagreen", lw=2, label=L"\sigma(\lambda_{2/3}, \tau)")
    plot!(0.001:0.001:10,-sch4, color="dodgerblue", lw=2, label=L"\sigma(\lambda_N, \tau)")
    plot!(0.001:0.001:10, schaefer_sigma, lw=4, color="navy", label=L"\sigma_{max}(\tau)")