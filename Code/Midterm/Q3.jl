using Random
using Plots

function sample_from_fair_die(rng)
    return rand(rng,1:6)
end


function sample_sum_of_rolls(rng)
    return sample_from_fair_die(rng) + sample_from_fair_die(rng)
end

function main()

    num_samples_array = collect(10:10:1000)
    prob_values = zeros(length(num_samples_array))  
    
    i = 1
    for num_samples in num_samples_array
        sum_array = zeros(12)
        for i in 1:num_samples
            rng = MersenneTwister()
            sampled_sum = sample_sum_of_rolls(rng)
            sum_array[sampled_sum] += 1
        end
        prob_values[i] = sum_array[8]/num_samples
        # println(sum_array)
        i+=1
    end

    
    snapshot = plot(size=(900,900), 
    dpi=300,
    # xticks=num_samples_array, 
    # yticks=0:0.1:0.4,
    xlabel="Number of Samples", 
    ylabel="Probability Value",
    title="Plot for Probability of getting the sum 8 with varying number of samples", 
    # axis=([], false),
    # legend=:bottom,
    legend=true,
    )

    true_prob = 5/36
    plot!(snapshot,num_samples_array,ones(length(num_samples_array))*true_prob,label="True Probability",linewidth=2)
    plot!(snapshot,num_samples_array,prob_values,label="Empirical Probability",linewidth=2)

    display(snapshot)
    return num_samples_array,prob_values,snapshot
end

#=
num_samples_array,prob_values,snapshot = main()
savefig(snapshot,"Q3_plot.png")
=#






