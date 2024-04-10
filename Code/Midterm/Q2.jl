given_data_matrix = 
[
    3.26 6.8 3.79 2.45 2.67 2.88 5.54 4.63 4.13 3.02 4.1 1.53 1.96 2.96 2.69 4.26 2.57 6.92 1.7 3.13 3.15 2.68 3.14 5.3 2.29 1.71 3.54 2.1 3.31 2.97 4.12 2.92 4.36 3.62 1.96 4.56 2.04 3.8 3.04 2.38 2.37 1.89 2.16 2.37 3.16 7.35 1.77 6.08 5.08 2.93 1.62 3.01 1.93 2.79 1.67 1.7 1.72 1.83 2.45 4.01 6.22 2.25 4.65 3.39 1.64 2.69 2.28 1.57 2.14 2.19 2.78 2.44 1.81 2.58 2.44 2.6 3.23 1.71 1.92 4.77 3.43 2.92 2.77 2.5 3.89 2.08 2.09 2.89 2.9 2.95 2.19 1.99 3.52 3.9 3.51 3.99 1.96 2.31 2.51 6.56 3 3.6 3.41 1.83 1.85 1.92 2.72 1.96 2.98 3.16 2.36 3.46 2.44 3.35 1.6 3.28 4.21 1.84 2.67 2.47 2.14 3.06 1.93 2.66 1.41 3.7 5.07 1.74 9.51 2.24 3.27 2.68 3.65 1.93 4.46 4.19 1.83 3.44 3.11 2.24 1.36 3.19 3.94 4.64 2.72 3.33 1.88 2.24 1.93 4.82 1.73 7.51 5.4 3.45 2.54 2.05 1.59 2.6 1.38 7.4 1.73 2.41 3.54 2.23 3.69 5.41 4.18 4.07 2.75 4.6 4.04 2.82 2.25 3 2.18 2.01 5.32 2.27 1.85 2.24 2.2 1.7 1.64 4.67 2.36 5.22 2.5 2.71 1.68 2.5 4.73 1.64 2.78 4.88 2.93 3.97 1.45 3.14 2.25 4.04
]

given_data = vec(given_data_matrix)
minimum_value = minimum(given_data)
start_value = floor(minimum_value)
maximum_value = maximum(given_data)
end_value = ceil(maximum_value)

num_bins = length(start_value:1.0:end_value)
counts = zeros(num_bins)

for value in given_data
    index = Int(floor(value))
    counts[index] += 1
end 

bin_values = collect(start_value:1.0:end_value)
bin_probability_values = counts/length(given_data)
dist = SparseCat(bin_values,bin_probability_values)
histogram(given_data, bins=num_bins, title="Histogram of Given Data", xlabel="Bins (intervals 1 unit apart)", ylabel="Frequency", legend=false)


num_samples = 1000
generated_data = []
for i in 1:num_samples
    sampled_bin = rand(dist)
    sampled_value = sampled_bin + rand()
    push!(generated_data,sampled_value)
end

generated_data_counts = zeros(num_bins)
for value in generated_data
    index = Int(floor(value))
    generated_data_counts[index] += 1
end
generated_data_bin_probability_values = generated_data_counts/num_samples
histogram(generated_data, bins=num_bins, title="Histogram of Generated Data", xlabel="Bins (intervals 1 unit apart)", ylabel="Frequency", legend=false)

# histogram(given_data, bins=num_bins, title="Histogram of given data", xlabel="Value", ylabel="Frequency", legend=false)