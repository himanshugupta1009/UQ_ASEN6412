include("spquad.jl")
include("utils.jl")
using StaticArrays



function do_smolyak_sparse_grid_clenshaw_curtis(num_dimensions, l, f)

    # Get the points and weights for the 1D quadrature
    points, weights = spquad(num_dimensions,l)
    points = [Tuple(c) for c in eachrow(points)]

    # Compute the integral
    start_index = 1
    integral = f(points[start_index])*weights[start_index]
    for i in (start_index+1):length(points)
        integral += (f(points[i])*weights[i])
    end

    return integral
end


function Q2_part3_mean(num_nodes,d,l)
    f(Y) = generate_u_sample(num_nodes,Y)
    integral = do_smolyak_sparse_grid_clenshaw_curtis(d, l, f)
    mean = 0.5*0.5*integral
    return mean
end

function compute_matrix_square(x)
    return x.*x
end

function Q2_part3_mean_square(num_nodes,d,l)
    f(Y) = compute_matrix_square(generate_u_sample(num_nodes,Y))
    integral = do_smolyak_sparse_grid_clenshaw_curtis(d, l, f)
    mean_square = 0.5*0.5*integral
    return mean_square
end


function visulization_Q2_part3()

    num_nodes = 100
    d = 2 #num_dimensions
    num_l = 4
    l_array = SVector(0,2,4,5)
    num_samples_array = MVector{num_l,Int}([2^i+1 for i in l_array])
    num_samples_array[1] = 1

    means = []
    variances = []

    for l in l_array
        M = Q2_part3_mean(num_nodes,d,l)
        M_square = Q2_part3_mean_square(num_nodes,d,l) 
        V = M_square .- compute_matrix_square(M)
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
