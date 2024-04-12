#=
By Counting the number of points on each axis, we can tell that there are 17 points.
This means l=4 because the number of points on each axis is 2^l+1.
So, we need to find 17 points in one dimension and take the tensor product of these points.
=#
include("spquad.jl")
using Plots


function get_grid_points(l1,l2)

    l1_one_d_points, l1_one_d_weights = spquad(1,l1)
    N1 = length(l1_one_d_points)  # N = 2^l + 1

    l2_one_d_points, l2_one_d_weights = spquad(1,l2)
    N2 = length(l2_one_d_points)  # N = 2^l + 1

    total_points = N1*N2
    grid_points = Array{NTuple{2,Float64},1}(undef,total_points)
    grid_weights = zeros(total_points)

    # println(N1, " ", N2, " ", total_points)

    index = 1
    for i in 1:N1
        for j in 1:N2
            # index = (i-1)*N1 + j
            # println(index, " ", i, " ", j)
            grid_points[index] = (l1_one_d_points[i],l2_one_d_points[j])
            grid_weights[index] = l1_one_d_weights[i]*l2_one_d_weights[j]
            index += 1
        end
    end

    return grid_points, grid_weights
end


function Q1_part1()

    q = 4
    d = 2

    # Possible values for l are 3 and 4
    possible_l1_l2_pairs_where_l_is_3 = SVector((3,0),(0,3),(1,2),(2,1)) 
    possible_l1_l2_pairs_where_l_is_4 = SVector((4,0),(0,4),(1,3),(3,1),(2,2))

    points_when_l_is_3 = Dict()
    for ele in possible_l1_l2_pairs_where_l_is_3
        l1, l2 = ele
        grid_points, grid_weights = get_grid_points(l1,l2)
        points_when_l_is_3[ele] = grid_points
    end
    points_when_l_is_4 = Dict()
    for ele in possible_l1_l2_pairs_where_l_is_4
        l1, l2 = ele
        grid_points, grid_weights = get_grid_points(l1,l2)
        points_when_l_is_4[ele] = grid_points
    end

    max_val = 1.1
    x_boundary = SVector{4,Float64}(-max_val,max_val,max_val,-max_val)
    y_boundary = SVector{4,Float64}(-max_val,-max_val,max_val,max_val)
    boundary = Shape(x_boundary, y_boundary)

    # Plot the points
    for ele in possible_l1_l2_pairs_where_l_is_3
        l1, l2 = ele
        grid_points = points_when_l_is_3[ele]
        total_points = length(grid_points)
        snapshot = plot(
            size=(700,700),
            dpi=300,
            xticks=-1.0:0.5:1.0, 
            yticks=-1.0:0.5:1.0,
            xlabel="y1", ylabel="y2",
            title="Tensor product of d=1 grids with l1=$l1 along x axis and \n 
                l2=$l2 along y axis with total points = $total_points \n",
            # axis=([], false),
            # legend=:bottom,
            legend=false
            # Distribution
            # Density
            )
        plot!(snapshot,boundary,linewidth=2,color=:white)
        for i in 1:total_points
            point = grid_points[i]
            scatter!(snapshot,point,linewidth=2,color=:green)
        end
        display(snapshot)
        savefig(snapshot, "./HW4/Q1_part1_l1=$l1"*"_l2=$l2.png")
    end

    for ele in possible_l1_l2_pairs_where_l_is_4
        l1, l2 = ele
        grid_points = points_when_l_is_4[ele]
        total_points = length(grid_points)
        snapshot = plot(
            size=(700,700),
            dpi=300,
            xticks=-1.0:0.5:1.0, 
            yticks=-1.0:0.5:1.0,
            xlabel="y1", ylabel="y2",
            title="Tensor product of d=1 grids with l1=$l1 along x axis and \n 
                l2=$l2 along y axis with total points = $total_points \n",
            # axis=([], false),
            # legend=:bottom,
            legend=false
            # Distribution
            # Density
            )
        plot!(snapshot,boundary,linewidth=2,color=:white)
        for i in 1:total_points
            point = grid_points[i]
            scatter!(snapshot,point,linewidth=2,color=:green)
        end
        display(snapshot)
        savefig(snapshot, "./HW4/Q1_part1_l1=$l1"*"_l2=$l2.png")
    end

end

#=
Q1_part1()
=#