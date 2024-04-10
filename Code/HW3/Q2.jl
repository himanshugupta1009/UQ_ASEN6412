import Roots as RT
import Distributions as DT
import Statistics as ST
using Random
using Plots

function get_λ(ω,l,σ)
    nume = 2*l
    deno = (l^2)*(ω^2) + 1
    λ = nume/deno
    return λ*σ*σ
end

#Need to find root of this equation
function odd_ω_equation(ω,a,l)
    return ω*tan(ω*a) - (1/l)
end 

#Need to find root of this equation
function even_ω_equation(ω,a,l)
    return ω + (1/l)*tan(ω*a)
end 

function odd_ϕ(t,ω,a)
    nume = cos(ω*t)
    deno = a
    deno += sin(2*ω*a)/(2*ω)
    deno = sqrt(deno)
    return nume/deno
end

function even_ϕ(t,ω,a)
    nume = sin(ω*t)
    deno = a
    deno -= sin(2*ω*a)/(2*ω)
    deno = sqrt(deno)
    return nume/deno
end

function tan_interval(i,a)
    # period of tan(ax) = pi/a
    δ = 0.00000001 # A really small number
    s = clamp( (i-1)*π/(2*a) + δ , 0.0, Inf)
    e = (i+1)*π/(2*a) - δ
    I = ( s,e )
    return I
end


function get_ωs(D,a,l)
    ωs = Array{Tuple{Int64,Float64},1}(undef,D)
    #=
    Stores both the index i and the value of ω.
    The index starts from 1 which is odd.
    Note, ω = 0.0 is a root of the even_ω_equation, but is not conisdred here because 
    when ω is 0.0, ϕ_even is also 0.0. So, it adds no value!
    =#
    for i in 1:D
        interval_id = 2*floor(i/2)
        I = tan_interval(interval_id,a)
        if( iseven(i) )
            ω = RT.find_zero(ω->even_ω_equation(ω,a,l),I)
            ωs[i] = (i,ω)
        else
            ω = RT.find_zero(ω->odd_ω_equation(ω,a,l),I)
            ωs[i] = (i,ω)
        end
    end
    return ωs
end


function generate_eigenpairs(D,a,l,σ)

    ωs = get_ωs(D,a,l)
    eigenpairs = Array{Tuple{Int64,Float64,Function},1}(undef,D)

    for i in 1:D
        (index,ω) = ωs[i]
        λ = get_λ(ω,l,σ)
        if( iseven(i) )
            eigenfunction_even(t) = even_ϕ(t,ω,a)
            eigenpairs[i] = (index,λ,eigenfunction_even)
        else
            eigenfunction_odd(t) = odd_ϕ(t,ω,a)
            eigenpairs[i] = (index,λ,eigenfunction_odd)
        end
    end
    return eigenpairs
end

#=
This function generates an entire trajectory of the GRP process.
Think of this as one trajectory with time for for a fixed ω.
=#
function generate_GRP_sample(D,a,l,mean,σ,max_T,Δt = 0.01,rng=MersenneTwister())
    
    #=
    Ideally, a shold be max_T/2.0
    =#

    eigenpairs = generate_eigenpairs(D,a,l,σ)
    Y = randn(rng,D)
    num_points = Int(max_T/Δt)+1
    X = Vector{Float64}(undef,num_points)
 
    for j in 1:num_points
        T = (j-1)*Δt - a
        s = 0.0 + mean
        for i in 1:D
            (index,λ,ϕ) = eigenpairs[i]
            s += sqrt(λ)*ϕ(T)*Y[i] 
        end
        X[j] = s
    end

    return X
end


function forcing_function(t)
    return -1.0
end

#Function that gives you Gd(x,ω) at any given x for a fixed ω determined by λ,ϕ and Y 
function KLE(t,eigenpairs,Y,mean,a)
    s = mean
    # Shift the time by a so that domain is [-a,a]
    T = t - a 
    for i in 1:length(eigenpairs)
        (index,λ,ϕ) = eigenpairs[i]
        s += sqrt(λ)*ϕ(T)*Y[i]
    end
    return s
end


#=
Solves Q2.Part b)
=#
include("../steady_state_math.jl")


function generate_u_sample(num_nodes,rng=MersenneTwister())

    D = 10
    a = 0.5
    l = 0.2
    mean = 1.0
    σ = 2
    eigenpairs = generate_eigenpairs(D,a,l,σ)
    Y = randn(rng,D)

    # num_nodes is n in the function fem1d_heat_steady
    slab_start = 0.0 #a
    slab_end = 1.0 #b
    boundary_condition_start = 0.0 #ua
    bc_end_distribution = DT.Normal(0.0,0.1) # N(μ=0.0, σ=0.1)
    boundary_condition_end = rand(rng,bc_end_distribution) #ub
    G(x) = KLE(x,eigenpairs,Y,mean,a)
    K(x) = exp(G(x)) #k
    F(x) = forcing_function(x) #f
    node_points = range(slab_start,stop=slab_end,length=num_nodes) #x

    sampled_u = fem1d_heat_steady(num_nodes,slab_start,slab_end,boundary_condition_start,
                                boundary_condition_end,K,F,node_points)

    return sampled_u
end


function do_monte_carlo_simulations(N,num_nodes,rng = MersenneTwister())
    s = rand(rng,UInt32)
    u_rng = MersenneTwister(s) 
    u_samples = [generate_u_sample(num_nodes) for i in 1:N]
    return u_samples
end


function main_Q2_part_b(N,num_nodes=100)
    u_samples = do_monte_carlo_simulations(N,num_nodes)
    M = ST.mean(u_samples)
    V = ST.var(u_samples)
    # println("Sample Mean of U samples over ",num_nodes," nodes : ",M)
    # println("Sample Variance of U samples: over ",num_nodes," nodes : ",V)
    return u_samples,M,V
end

function visulization_Q2_partb()

    num_nodes = 100
    num_samples_array = (100,200,500,1000,2000,5000,10000,20000,50000)
    # num_samples_array = (100,200)
    means = []
    variances = []
    for N in num_samples_array
        u_samples,M,V = main_Q2_part_b(N,num_nodes)
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
s1,s2 = visulization_Q2_partb()
Plots.savefig(s1,"HW3/q2_mean_convergence.png")
Plots.savefig(s2,"HW3/q2_variance_convergence.png")
=#

#=
Solves Q2.Part c)
=#
function find_u_max_mean_covar(u_samples)
    u_max = [maximum(u) for u in u_samples]
    mean_u_max = ST.mean(u_max)
    var_u_max = ST.var(u_max)
    return mean_u_max,var_u_max
end


function main_Q2_part_c()
    N = 10000
    num_nodes = 100
    u_samples = do_monte_carlo_simulations(N,num_nodes)
    mean_u_max,var_u_max = find_u_max_mean_covar(u_samples)

    num_successful_samples = 0
    threshold = mean_u_max + 3*sqrt(var_u_max)
    u_samples = do_monte_carlo_simulations(N,num_nodes)
    for u in u_samples
        u_max = maximum(u)
        if(u_max > threshold)
            num_successful_samples += 1
        end
    end

    println("Number of samples that satisfied the condition : ",num_successful_samples)
    println("Desired Probability : ",num_successful_samples/N)
    return num_successful_samples

end
#=
main_Q2_part_c()    
=#