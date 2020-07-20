using CSV, DataFrames

## Read in data & tree structure

data_path = raw"C:\Users\gabeh\OneDrive\Documents\Biomath\Julia tutorial\data\yeast_continuous.csv"
tree_path = raw"C:\Users\gabeh\OneDrive\Documents\Biomath\Julia tutorial\data\yeast.txt"

data_df = CSV.read(data_path)
taxa = String.(data_df.taxon)
traits = Matrix(data_df[!, 2:end]) # first column stores the taxa, the rest store the actual trait data
                                   # I should really check to make sture that the CSV file is properly structured

tree = read(tree_path, String)

## Set model parameters
k = 2
shape_mult = 10.0^4
run_name = "old_test"

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


xml_path = joinpath(@__DIR__, "$run_name.xml")
make_xml(xml_path, taxa, traits, tree, k, shape_mult)
