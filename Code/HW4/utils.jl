include("../steady_state_math.jl")
include("../helper_functions.jl")

#Function that gives you K(x,ω) at any given x for a fixed ω determined by λ,ϕ and Y 
function calculate_K_xω(x,eigenpairs,Y)
    K_bar = 1
    σ_K = 0.1
    d = length(eigenpairs)
    a = 0.5   #IS this needed? Should it be zero for this Q?

    s = 0.0
    T = x - a
    for i in 1:length(eigenpairs)
        (index,λ,ϕ) = eigenpairs[i]
        s += sqrt(λ)*ϕ(T)*Y[i]
    end
    s = K_bar + σ_K*s
    return s
end

function forcing_function(t)
    return -1.0
end

function generate_u_sample(num_nodes,Y)

    D = 2
    #D = length(Y)
    a = 0.5
    l = 2.0
    σ = 1
    eigenpairs = generate_eigenpairs(D,a,l,σ)

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