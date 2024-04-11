#=
By Counting the number of points on each axis, we can tell that there are 17 points.
This means l=4 because the number of points on each axis is 2^l+1.
So, we need to find 17 points in one dimension and take the tensor product of these points.
=#
include("spquad.jl")
using Plots

function Q1_part1()

    l = 4
    d = 2

    one_d_points, one_d_weights = spquad(1,l)
    N = length(one_d_points)  # N = 2^l + 1
    total_points = N^d

    # println("Number of points on each axis: ", length(one_d_points))
    # println("Number of points in total: ", length(one_d_points)^2)

    grid_points = Array{NTuple{d,Float64},1}(undef,total_points)
    grid_weights = zeros(total_points)

    for i in 1:N
        for j in 1:N
            index = (i-1)*N + j
            grid_points[index] = (one_d_points[i],one_d_points[j])
            grid_weights[index] = one_d_weights[i]*one_d_weights[j]
        end
    end

    snapshot = plot(size=(700,700), dpi=300,
        # xticks=:0.1:1.0, yticks=0:0.2:3,
        xlabel="y1", ylabel="y2",
        title="Tensor product of d=1 grids with l=4 \n (i.e. 17 points) for each axis",
    # axis=([], false),
    # legend=:bottom,
    legend=false
    # Distribution
    # Density
    )

    for i in 1:total_points
        point = grid_points[i]
        scatter!(snapshot,point,linewidth=2)
    end
    display(snapshot)
    return snapshot
end

#=
Q1_part1()
=#