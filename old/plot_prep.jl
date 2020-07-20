using CSV, DataFrames, Statistics
using BeastUtils.Logs

function process_log(log_path::String, csv_path::String, data_path::String,
                     labels_path::String; burnin = 0.1, keep_threshold = 0.9)

    cols, data = get_log(log_path, burnin= burnin)
    L_header = "L"
    sv_header = "sv"

    L_inds = findall(x -> startswith(x, L_header), cols)
    sv_inds = findall(x -> startswith(x, sv_header), cols)
    L_cols = cols[L_inds]
    L_data = @view data[:, L_inds]

    k = length(sv_inds)
    p, r = divrem(length(L_cols), k)
    @assert r == 0

    n = size(data, 1)
    upper_threshold = Int(floor(keep_threshold * n))
    lower_threshold = Int(ceil((1.0 - keep_threshold) * n))

    L = Matrix{Union{Missing, Float64}}(undef, k, p)
    fill!(L, missing)

    labels_df = CSV.read(labels_path)
    labels = labels_df.label
    new_names = labels_df.pretty
    trait_types = labels_df.cat

    original_labels = string.(names(CSV.read(data_path))[2:end]) # TODO: this is super inefficient, just get labels

    @assert length(labels) == length(original_labels)
    perm = indexin(original_labels, labels)

    @assert length(findall(isnothing, perm)) == 0

    for i = 1:k
        for j in 1:p
            col = (i - 1) * p + j
            @assert L_cols[col] == "$L_header$i$(j)"
            vals = @view(L_data[:, col])
            n_pos = count(x -> x > 0.0, vals)
            if n_pos > upper_threshold || n_pos < lower_threshold
                L[i, perm[j]] = mean(@view L_data[:, col])
            end
        end
    end

    df = DataFrame()
    df.L = vec(L)
    df.col = repeat(1:p, inner=k)
    df.row = repeat(1:k, outer=p)

    trait_names = string.(names(CSV.read(data_path))[2:end])
    df.trait = repeat(new_names, inner=k)
    df.cat = repeat(trait_types, inner=k)
    CSV.write(csv_path, df)
end

labels_path = raw"C:\Users\gabeh\OneDrive\Documents\Biomath\Julia tutorial\data\yeast_labels.csv"
labels_df = CSV.read(labels_path)
cat_levs = unique(labels_df.cat)

csv_path = "old_test.csv"
process_log(processed_path, csv_path, data_path, labels_path)

plot_path = "old_test.pdf"
fact = 1:k

# below needed to avoid issues with '°' character for temperatures
tmp_path = "tmp.csv"
CSV.write(tmp_path, DataFrame(levs = labels_df.pretty))
