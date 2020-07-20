using Blink, Interact # GUI application
cd(@__DIR__) # change working directory to directory of this file

## GUI


name_field = textbox(label="Name this BEAST run: ")
data_import = filepicker("Import Data") # button for picking the data file
tree_import = filepicker("Import Tree") # button for picking the tree file
labels_import = filepicker("Import Plotting Labels") # button for picking file with plot labels
k_box = spinbox(1:10, label="Number of factors: ") # input for selecting the number of factors
shrink_box = spinbox(1:10, label="Shrinkage intensity: ") # input for selecting the appropriate shrinkage
start_button = button("Start Run") # button for starting run

data_alert = alert("Please select an existing data file before running.")
tree_alert = alert("Please select an existing tree file before running.")
labels_alert = alert("Please select an existing labels file before running.")

done = false
w = Window()
body!(w, vbox(name_field, data_import, tree_import, labels_import, k_box, shrink_box, start_button, data_alert, tree_alert))

on(start_button) do _ # what to do when the `start_button` is pressed
    if !isfile(data_import[]) # check that a file for the data has been provided
        data_alert()
    elseif !isfile(tree_import[]) # check that a file for the tree has been provided
        tree_alert()
    elseif !isfile(labels_import[]) # check that a file for the tree has been provided
        labels_alert()
    else
        global done = true
    end
end

while !done # this is a really bad way to do this, but I'm not good enough to find another way
    if !isopen(w.content.sock)
        error("Closed window without running analysis. Re-run this script to start over.")
    end
    sleep(1.0)
end

close(w)


## Read in data & tree structure

using CSV, DataFrames

data_path = data_import[]
tree_path = tree_import[]
labels_path = labels_import[]

data_df = CSV.read(data_path)
taxa = String.(data_df.taxon)
traits = Matrix(data_df[!, 2:end]) # first column stores the taxa, the rest store the actual trait data
                                   # I should really check to make sture that the CSV file is properly structured

tree = read(tree_path, String)

## Set model parameters
k = k_box[]
shape_mult = 10.0^(shrink_box[])
run_name = name_field[]

## Build XML file (the instructions file for BEAST)

using BeastUtils.XMLConstructor

function make_xml(path::String,
                  taxa::Vector{String},
                  data::Matrix{Float64},
                  tree::String,
                  k::Int,
                  shape_mult::Float64;
                  chain_length::Int = 10000)

    shapes = fill(shape_mult, k - 1)

    bx = XMLConstructor.make_PFA_XML(data, taxa, tree, k,
                                     useHMC = false,
                                     shrink_loadings = true,
                                     chain_length = chain_length,
                                     fle = div(chain_length, 1000),
                                     sle = 1000)

    facs = XMLConstructor.get_integratedFactorModel(bx)
    XMLConstructor.set_shrinkage_mults!(facs, shapes = shapes)
    XMLConstructor.save_xml(path, bx)
    return nothing
end

if !isdir(run_name)
    mkdir(run_name)
end

cd(run_name)

xml_path = "$run_name.xml"
make_xml(xml_path, taxa, traits, tree, k, shape_mult)

## Run BEAST from Julia using the Cmd type

beast_path = joinpath(@__DIR__, "beast.jar")

command = Cmd(["java", "-jar", beast_path, "-overwrite", xml_path]) # same as typing `java -jar <path_to_beast_jar_file> -overwrite <my_xml_file>`
                                                                    # into the command line
run(command)


## Data post-processing

using SVDProcessing

log_path = "$run_name.log"
processed_path = "$(run_name)_svd.log"
SVDProcessing.svd_logs(log_path, processed_path)

## Plotting

using RCall, Statistics
using BeastUtils.Logs


R"""
library(ggplot2)
plot_loadings <- function(csv_path, plot_name, trait_levs, cat_levs){
    df  <- read.csv(csv_path, header=TRUE)

    df$trait <- factor(df$trait, levels=trait_levs)
    df$cat <- factor(df$cat, levels=cat_levs)
    df$L <- sapply(df$L, as.numeric)
    ggplot(df, aes(x=trait, y=row, fill=L)) +
            facet_grid(~ cat, scales="free_x", space="free_x") +
            geom_tile() +
            scale_fill_gradient2(low="orange", mid="white", high="purple", midpoint=0) +
            scale_x_discrete(position = "top") +
            scale_y_discrete(limits=fact) +
      labs(y="Loadings Row", x="Trait") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle=90, hjust=0),
            strip.text.x = element_blank()
            )
            # axis.title.y = element_text())

    ggsave(plot_name, width=15, height=10, units="in")
    gc()
}

plot_factors <- function(csv_path, plot_name){
    df  <- read.csv(csv_path, header=TRUE)

    ggplot(df, aes(x=f1, y=f2, color=class)) +
            geom_point() +
            labs(y="f2", x="f1") +
      theme_minimal()

    ggsave(plot_name, width=6, height=4, units="in")

gc()
}
"""

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


labels_df = CSV.read(labels_path)
cat_levs = unique(labels_df.cat)

csv_path = "$run_name.csv"
process_log(processed_path, csv_path, data_path, labels_path)

@rput csv_path
plot_path = "$run_name.pdf"
@rput plot_path
@rput cat_levs
fact = 1:k
@rput fact

# below needed to avoid issues with '°' character for temperatures
tmp_path = "tmp.csv"
@rput tmp_path
CSV.write(tmp_path, DataFrame(levs = labels_df.pretty))

R"""
pretty_names <- read.csv(tmp_path)$levs
plot_loadings(csv_path, plot_path, pretty_names, cat_levs)
"""

rm(tmp_path)
