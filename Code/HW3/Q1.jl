using Random
import Statistics as ST
import Roots as RT 

function inverse_F1(t)
    #=
    F1(x) = 1-e^(-x)
    t = 1-e^(-x)
    x = -ln(1-t)
    =#
    return -log(1-t)
end

function sample_X1(rng=MersenneTwister())
    z1 = rand(rng)
    x1 = inverse_F1(z1)
    return x1
end

function F2_given_x1(x1,x2)
    value = 1 - (1+x2)*exp(-(1+x1)*x2)
    return value
end

function inverse_F2_given_x1(t,x1)
    #=
    F2(x2|x1) = 1 - (1+x2)*e^(-(1+x1)*x2)
    t = F2(x2|x1) 
    0 = t - F2(x2|x1)
    0 = 1 - t - (1+x2)*e^(-(1+x1)*x2)
    Find roots of this equation, i.e. Find x2 such that t - F2(x2|x1) = 0
    Choose the positive one!
    =#
    g(x) = t - F2_given_x1(x1,x)
    interval = (0.0,Inf)
    root = RT.find_zero(g,interval)
    return root
end

function sample_X2_given_x1(x1,rng=MersenneTwister())
    z2 = rand(rng)
    x2 = inverse_F2_given_x1(z2,x1)
    return x2
end

function generate_samples_X1X2(N)
    rng = MersenneTwister()
    sampled_X1 = [sample_X1(rng) for i in 1:N]
    sampled_X2 = sample_X2_given_x1.(sampled_X1,rng)
    return sampled_X1,sampled_X2
end


function generate_samples_kzkα(N)
    sampled_X1,sampled_X2 = generate_samples_X1X2(N)
    kz0 = 1
    kα0 = 10
    sampled_kz = kz0 .+ sampled_X1
    sampled_kα = kα0 .+ sampled_X2
    return sampled_kz,sampled_kα
end

function main_Q1(N=100000)
    sampled_X1,sampled_X2 = generate_samples_X1X2(N)
    sampled_X1X2 = sampled_X1.*sampled_X2
    println("Sample Mean of X1: ",ST.mean(sampled_X1))
    println("Sample Mean of X2: ",ST.mean(sampled_X2))
    println("Sample Mean of X1*X2: ",ST.mean(sampled_X1X2))
    return sampled_X1,sampled_X2
end

main_Q1();