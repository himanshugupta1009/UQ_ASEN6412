import Roots as RT
import Distributions as DT
import Statistics as ST
using Random

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


function find_D(a,l,σ,α=0.9,k=4)

    D0 = 1
    while(true)
        D = D0+k
        eigenpairs = generate_eigenpairs(D,a,l,σ)
        numerator = 0.0
        denominator = 0.0
        for i in 1:D
            (index,λ,ϕ) = eigenpairs[i]
            if(i <= D0)
                numerator += λ
                denominator += λ
            else
                denominator += λ
            end 
        end
        if(numerator/denominator >= α)
            println("Numerator Value is : $numerator")
            println("Denominator Value is : $denominator")
            println("The Ratio is : ",numerator/denominator)
            return D0
        else
            D0 += 1
        end
    end

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


#=
This function also generates an entire trajectory of the GRP process.
But it assumes the eigenpairs are given.
Think of this as one trajectory with time for for a fixed ω.
=#
function generate_GRP_sample(eigenpairs,mean,a,max_T,Δt = 0.01,rng=MersenneTwister())
    
    D = length(eigenpairs)
    Y = randn(rng,D)
    num_points = Int(max_T/Δt)+1
    # X = Vector{Float64}(undef,num_points)
    X = Matrix{Float64}(undef,num_points,1)
 
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

function main_Q2(l)
    
    #Define all the necessary variables
    a = 0.5
    σ = 2.0
    α = 0.9
    k = 1000 # Value for checking if the threshold α is met
    mean = 1.0
    max_T = 1.0  #For this problem, this variable is actually max_X
    Δt = 0.01    #For this problem, this variable is actually Δx


    #Part A)
    D = find_D(a,l,σ,α,k)
    println("The value of D is : $D")
    eigenpairs = generate_eigenpairs(D,a,l,σ)

    #Part B)
    num_samples = 3
    G_samples = [generate_GRP_sample(eigenpairs,mean,a,max_T,Δt) for i in 1:num_samples]
    snapshot1 = plot(size=(700,700), dpi=300,
    # xticks=:0.1:1.0, yticks=0:0.2:3,
    xlabel="x", ylabel="G(x,ω)",
    title="Realizations of G(x,ω)" 
    # axis=([], false),
    # legend=:bottom,
    # legend=false
    # Distribution
    # Density
    )
    for i in 1:num_samples
        plot!(snapshot1,0:Δt:max_T,G_samples[i],label="ω$i",linewidth=2)
    end
    display(snapshot1)

    #Part C)
    num_samples = 100000
    G_samples = [generate_GRP_sample(eigenpairs,mean,a,max_T,Δt) for i in 1:num_samples]
    mean_G_samples = ST.mean(G_samples)
    variance_G_samples = ST.var(G_samples)
    
    snapshot_mean = plot(size=(700,700), dpi=300,
    # xticks=:0.1:1.0, yticks=0:0.2:3,
    xlabel="x", ylabel="G(x,ω)",
    title="Mean value of G(x,ω) across $num_samples samples" 
    # axis=([], false),
    # legend=:bottom,
    # legend=false
    # Distribution
    # Density
    )
    plot!(snapshot_mean,0:Δt:max_T,mean_G_samples,label="Mean G(x,ω)",linewidth=2)
    plot!(snapshot_mean,0:Δt:max_T,[mean for i in 0:Δt:max_T],label="True Mean G(x,ω)",linewidth=2)
    display(snapshot_mean)

    snapshot_var = plot(size=(700,700), dpi=300,
    # xticks=:0.1:1.0, yticks=0:0.2:3,
    xlabel="x", ylabel="G(x,ω)",
    title="Variance of G(x,ω) across $num_samples samples" 
    # axis=([], false),
    # legend=:bottom,
    # legend=false
    # Distribution
    # Density
    )
    plot!(snapshot_var,0:Δt:max_T,variance_G_samples,label="Variance G(x,ω)",linewidth=2)
    plot!(snapshot_var,0:Δt:max_T,[σ*σ for i in 0:Δt:max_T],label="True Variance G(x,ω)",linewidth=2)
    display(snapshot_var)

    return snapshot1,snapshot_mean,snapshot_var
end

#Q3 - Part 1
l1 = 2.0
s1,s2,s3 = main_Q2(l1)
savefig(s1,"HW2/q3_part1_G_samples.png")
savefig(s2,"HW2/q3_part1_mean_G_samples.png")
savefig(s3,"HW2/q3_part1_variance_G_samples.png")

#Q3 - Part 2
l2 = 0.2
s1,s2,s3 = main_Q2(l2)
savefig(s1,"HW2/q3_part2_G_samples.png")
savefig(s2,"HW2/q3_part2_mean_G_samples.png")
savefig(s3,"HW2/q3_part2_variance_G_samples.png")
