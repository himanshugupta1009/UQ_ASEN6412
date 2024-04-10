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

