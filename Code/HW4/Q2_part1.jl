include("utils.jl")
using Plots
using StaticArrays
import Distributions as DT

#=
Q2 part 1
=#
function sample_LHS_points(num_partitions,num_dimensions,seeds_array)

    #=
        num_partitions is essentially the value of N in each dimension.
    =#
    
    dim_start = -1.0
    dim_end = 1.0  
    Δx = (dim_end-dim_start)/(num_partitions)
    dim_sampled_points = zeros(num_partitions)

    sampled_points = Dict{Int,Array{Float64,1}}()

    for d in 1:num_dimensions
        rng = MersenneTwister(seeds_array[d])
        dim_sampled_points = zeros(num_partitions)
        for i in 1:num_partitions
            a = dim_start + (i-1)*Δx
            b = dim_start + i*Δx
            point = rand(rng, DT.Uniform(a,b))
            dim_sampled_points[i] = point
            # println((a,b))
        end
        sampled_points[d] = dim_sampled_points
    end

    return sampled_points
end
#=
seeds_array = (11,111)
num_partitions = 10
num_dimensions = 2
sample_LHS_points(num_partitions,num_dimensions,seeds_array)
=#


function generate_Ys(num_samples,num_dimensions,rng=MersenneTwister(1))

    num_partitions = num_samples
    seeds_array = (11,111)
    sampled_points = sample_LHS_points(num_partitions,num_dimensions,seeds_array)

    Ys = Array{NTuple{num_dimensions,Float64},1}(undef,num_samples)
    shuffled_indices = Dict{Int,Array{Int64,1}}()

    for i in 1:num_dimensions
        shuffled_indices[i] = shuffle(rng,1:num_partitions)
    end

    for j in 1:num_samples
        sampled_Y = MVector{num_dimensions,Float64}(undef)
        for d in 1:num_dimensions
            sampled_Y[d] = sampled_points[d][shuffled_indices[d][j]]
        end
        Ys[j] = Tuple(sampled_Y)
    end

    return Ys
end


function do_LHS_sampling(N,num_dimensions,num_nodes,rng = MersenneTwister(9))
    Ys = generate_Ys(N,num_dimensions,rng)
    u_samples = [generate_u_sample(num_nodes,Ys[i]) for i in 1:N]
    return u_samples
end


function Q2_part_1(N,num_nodes=100,rng = MersenneTwister(19))
    num_dimensions = 2
    u_samples = do_LHS_sampling(N,num_dimensions,num_nodes,rng)
    M = ST.mean(u_samples)
    V = ST.var(u_samples)
    # println("Sample Mean of U samples over ",num_nodes," nodes : ",M)
    # println("Sample Variance of U samples: over ",num_nodes," nodes : ",V)
    return u_samples,M,V
end


function visulization_Q2_part1()

    num_nodes = 100
    num_dimensions = 2
    num_samples_array = (100,200,500,1000,2000,5000)
    # num_samples_array = (100,200,500,1000,2000,5000,10000,20000,50000)
    # num_samples_array = (100,200)
    # num_samples_array = (1:1:10)
    means = []
    variances = []
    seed = 13
    for N in num_samples_array
        rng = MersenneTwister(seed)
        u_samples,M,V = Q2_part_1(N,num_nodes,rng)
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

#=
num_samples = 10
num_nodes = 101
Q2_part_1(num_samples,num_nodes,MersenneTwister(11))
visulization_Q2_part1()
=#