import StatsBase as SB
using Plots
import KernelDensitySJ as KD

function pdf_chebyshev(x)
    v = (1-x^2)
    v = sqrt(1/v) * (1/pi)
    return v
end

function cdf_chebyshev(x)
    v = (1/pi)*asin(x)
    v += 0.5
    return v
end

function sample_x(z)
    v = pi*(z-0.5)
    return sin(v)
end

function generate_samples_X(N)

    sampled_U = rand(N)
    sampled_X = sample_x.(sampled_U)

    return sampled_U,sampled_X
end

function main_Q2()

    N = 10000
    sampled_U,sampled_X = generate_samples_X(N)
    
    X_points = -1:0.001:1

    #Empirical Density
    true_pdf_values = pdf_chebyshev.(X_points)
    bw = KD.bwsj(sampled_X)
    empirical_pdf_values = KD.density(sampled_X,bw,X_points)
    p = visualize(X_points, true_pdf_values, "True PDF Values")
    plot!(X_points,empirical_pdf_values,linewidth=2,color=:green,label="Emprical PDF Values")

    #Distribution function
    true_cdf_values = cdf_chebyshev.(X_points)
    e = SB.ecdf(sampled_X)
    empirical_CDF_values = e.(X_points)
    p = visualize(X_points, true_cdf_values, "True CDF Values")
    plot!(X_points,empirical_CDF_values,linewidth=2,color=:green,label="Emprical CDF Values")    
    return
end

function visualize(x_data, y_data, lab)

    snapshot = plot(size=(700,700), dpi=300,
        xticks=-1.0:0.1:1.0, yticks=0:0.2:3,
        xlabel="X", ylabel="Probability Distribution Value",
        title="Probability Distribution Function" 
        # axis=([], false),
        # legend=:bottom,
        # legend=false
        # Distribution
        # Density
        )
    x = collect(-1.0:0.1:1.0)
    plot!(snapshot,x_data,y_data,linewidth=2,color=:red,label=lab)
    return snapshot
end
