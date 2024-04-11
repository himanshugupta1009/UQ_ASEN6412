include("Q2_part1.jl")
include("Q2_part2.jl")
include("Q2_part3.jl")

#=
Do MC
Do LHS
Do Tensor Product
Do Smolyak with CC
=#


function generate_u_sample_MC(num_nodes,rng=MersenneTwister(7))

    D = 2
    a = 0.5
    l = 2.0
    σ = 1
    eigenpairs = generate_eigenpairs(D,a,l,σ)
    Y = randn(rng,D)

    # num_nodes is n in the function fem1d_heat_steady
    slab_start = 0.0 #a
    slab_end = 1.0 #b
    boundary_condition_start = 0.0 #ua
    boundary_condition_end = 0.0 #ub
    K(x) = calculate_K_xω(x,eigenpairs,Y) #k
    F(x) = forcing_function(x) #f
    node_points = range(slab_start,stop=slab_end,length=num_nodes) #x

    sampled_u = fem1d_heat_steady(num_nodes,slab_start,slab_end,boundary_condition_start,
                                boundary_condition_end,K,F,node_points)

    return sampled_u
end

function do_monte_carlo_simulations(N,num_nodes,rng = MersenneTwister(77))
    u_samples = [generate_u_sample_MC(num_nodes,rng) for i in 1:N]
    return u_samples
end

function main_MC(N,num_nodes=100,rng=MersenneTwister(777))
    u_samples = do_monte_carlo_simulations(N,num_nodes,rng)
    M = ST.mean(u_samples)
    V = ST.var(u_samples)
    # println("Sample Mean of U samples over ",num_nodes," nodes : ",M)
    # println("Sample Variance of U samples: over ",num_nodes," nodes : ",V)
    return u_samples,M,V
end

function Q2_part4()

    num_dimensions = 2
    num_nodes = 101
    x_index = 51 #index of x = 0.5
    #For MC and LHS
    num_samples_array = collect(0:1:128) 
    num_samples_array[1] = 1
    push!(num_samples_array,129)
    #For TP-CC and Smolyak-CC
    l_array = SVector(0,1,2,3,4,5,6,7)
    num_l = length(l_array)
    num_samples_array_using_l = MVector{num_l,Int}([2^i+1 for i in l_array])
    num_samples_array_using_l[1] = 1


    # MC
    means_MC = []
    seed = rand(UInt32)
    for N in num_samples_array
        rng = MersenneTwister(seed)
        u_samples,M,V = main_MC(N,num_nodes,rng)
        push!(means_MC,M[51])
    end

    #LHS
    means_LHS = []
    seed = 13
    for N in num_samples_array
        rng = MersenneTwister(seed)
        u_samples,M,V = Q2_part_1(N,num_nodes,rng)
        push!(means_LHS,M[51])
    end

    #Tensor Product
    means_TP = []
    d = num_dimensions
    for l in l_array
        M = Q2_part2_mean(num_nodes,d,l)
        push!(means_TP,M[51])
    end

    #Smolyak
    means_SG = []
    d = num_dimensions
    for l in l_array
        M = Q2_part3_mean(num_nodes,d,l)
        push!(means_SG,M[51])
    end

    return means_MC,means_LHS,means_TP,means_SG,num_samples_array,
                    num_samples_array_using_l

end


function visulization_Q2_part4()


    means_MC,means_LHS,means_TP,means_SG,num_samples_array,
    num_samples_array_using_l = Q2_part4()


    snapshot = plot(size=(700,700), dpi=300,
        # xticks=:0.1:1.0, yticks=0:0.2:3,
        xlabel="Number of Samples", ylabel="Mean of 'u' value at x=0.5",
        title="Plot for Mean of 'u' at x=0.5 across different methods with \n varying number of samples" 
    # axis=([], false),
    # legend=:bottom,
    # legend=false
    # Distribution
    # Density
    )

    plot!(snapshot,num_samples_array,means_MC,label="Monte Carlo",linewidth=2)
    plot!(snapshot,num_samples_array,means_MC,label="LHS",linewidth=2)
    plot!(snapshot,num_samples_array_using_l,means_TP,label="Tensor Product",linewidth=2)
    plot!(snapshot,num_samples_array_using_l,means_TP,label="Smolyak-CC",linewidth=2)

    display(snapshot)
    return snapshot

end

#=
Q2_part4()
visulization_Q2_part4()
=#


function visulization_just_MC()

    num_nodes = 101
    # num_samples_array = (100,200,500,1000,2000,5000,10000,20000,50000)
    num_samples_array = (1:1:5)
    means = []
    variances = []
    for N in num_samples_array
        u_samples,M,V = main_MC(N,num_nodes)
        push!(means,M)
        push!(variances,V)
    end

    snapshot1 = plot(size=(700,700), dpi=300,
        # xticks=:0.1:1.0, yticks=0:0.2:3,
        xlabel="Node i", ylabel="Mean of 'u' value for Node i ",
        title="Plot for Mean of 'u' value across different nodes with \n varying number of samples" 
    # axis=([], false),
    # legend=:bottom,
    # legend=false
    # Distribution
    # Density
    )

    for i in 1:length(num_samples_array)
        M = means[i]
        N = num_samples_array[i]
        plot!(snapshot1,1:num_nodes,M,label="N = $N",linewidth=2)
    end
    display(snapshot1)

    snapshot2 = plot(size=(700,700), dpi=300,
    # xticks=:0.1:1.0, yticks=0:0.2:3,
    xlabel="Node i", ylabel="Variance of 'u' value for Node i ",
    title="Plot for Variance of 'u' value across different nodes with \n varying number of samples" 
    # axis=([], false),
    # legend=:bottom,
    # legend=false
    # Distribution
    # Density
    )

    for i in 1:length(num_samples_array)
        V = variances[i]
        N = num_samples_array[i]
        plot!(snapshot2,1:num_nodes,V,label="N = $N",linewidth=2)
    end
    display(snapshot2)

    return snapshot1,snapshot2
end