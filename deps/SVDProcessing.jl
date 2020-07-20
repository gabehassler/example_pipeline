module SVDProcessing

using LinearAlgebra, Statistics, LinearAlgebra.BLAS
using BeastUtils.Logs, BeastUtils.MatrixUtils, BeastUtils.ESS


const DEBUG = false

function process_log(log_path::String, header::Union{Regex, String}, k::Int,
                    p::Int;
                    fac_header::String = "factors.",
                    rotate_factors::Bool = false,
                    find_best::Bool = true,
                    rotation_cols::Vector{Int} = zeros(Int, k),
                    transposed::Bool = false)

    cols, data = get_log(log_path, burnin = 0.0)
    L_inds = findall(x -> startswith(x, header), cols)
    @assert length(L_inds) == k * p

    f_inds = rotate_factors ? findall(x -> startswith(x, fac_header), cols) : Int[]
    n_taxa = div(length(f_inds), k)
    @assert length(f_inds) == n_taxa * k

    L_cols = cols[L_inds]
    L_data = data[:, L_inds]

    f_cols = cols[f_inds]
    f_data = data[:, f_inds]



    # rotation_cols = ones(Int, k)
    # if find_best
    #     rotation_cols .= find_rotation_cols(data, k, p)
    # end


    # data = randn(10000, k * p)
    # cols = ["col$i" for i = 1:(k * p)]

    n = size(data, 1)
    # n = 100 #TODO: remove


    d_storage = zeros(n, k)
    v_storage = zeros(n, k * p)
    l_storage = zeros(n, k * p)
    f_storage = zeros(n, n_taxa * k)


    L = zeros(k, p)
    Ft = zeros(k, n_taxa)
    UtFt = zeros(k, n_taxa)

    for i = 1:n
        fill_L!(L, L_data, i, transposed)

        s = svd(L)
        for j = 2:k
            @assert s.S[j - 1] >= s.S[j]
        end
        # rotate_svd!(s, cols = rotation_cols)
        if DEBUG
            L2 = s.U * Diagonal(s.S) * s.Vt
            @assert maximum(abs.(L2 - L)) < 1e-12
        end
        d_storage[i, :] .= s.S
        v_storage[i, :] .= vec(s.Vt')

        l_storage[i, :] .= vec(s.Vt' * Diagonal(s.S))

        if rotate_factors
            copyto!(Ft, @view f_data[i, :])
            gemm!('T', 'N', 1.0, s.U, Ft, 0.0, UtFt)
            f_storage[i, :] .= vec(UtFt)
        end
    end

    force_rotation = false
    if any(rotation_cols .!= 0)
        force_rotation = true
    end
    if !force_rotation
        if !find_best
            fill!(rotation_cols, 1)
        else
            rotation_cols .= find_rotation_cols(v_storage, k, p)
        end
    end
    reflect!(v_storage, l_storage, f_storage, rotation_cols, k, p)
    return d_storage, v_storage, l_storage, f_storage
end

function fill_L!(L::Matrix{Float64}, data::Matrix{Float64}, row::Int, transposed::Bool)
    k, p = size(L)
    if transposed
        for i = 1:p
            offset = k * (i - 1)
            for j = 1:k
                L[j, i] = data[row, offset + j]
            end
        end
    else
        for i = 1:k
            offset = p * (i - 1)
            for j = 1:p
                L[i, j] = data[row, offset + j]
            end
        end
    end
end

function rotate_svd!(s::SVD{Float64,Float64,Array{Float64,2}}; cols::Vector{Int} = ones(Int, length(s.S)))
    k = length(s.S)
    rev = ones(k)
    for i = 1:k
        if s.Vt[i, cols[i]] < 0.0
            rev[i] = -1.0
        end
    end
    Q = Diagonal(rev)
    s.U .= s.U * Q #TODO: make memory efficient
    s.Vt .= Q * s.Vt #TODO: see above
end

function find_rotation_cols(data::Matrix{Float64}, k::Int, p::Int)
    abs_mat = abs.(data)
    means = vec(mean(abs_mat, dims = 1))
    vars = vec(var(abs_mat, dims = 1))

    zs = means ./ sqrt.(vars)
    nan_inds = findnans(zs)
    zs[nan_inds] .= 0.0
    cols = zeros(Int, k)
    for i = 1:k
        l = ((i - 1) * p) + 1
        u = i * p
        zmax = findmax(zs[l:u])
        cols[i] = zmax[2]
    end
    return cols
end

function process_loadings(log_path::String, header::Union{Regex, String}, k::Int, p::Int; find_best::Bool = true)
    cols, data = get_log(log_path, header, burnin = 0.0)

    @assert length(cols) == k * p

    rotation_cols = ones(Int, k)
    if find_best
        rotation_cols .= find_rotation_cols(data, k, p)
    end
    reflect!(data, rotation_cols, k, p)
    return data, cols
end

function reflect!(data::Matrix{Float64}, loadings::Matrix{Float64},
                    factors::Matrix{Float64}, cols::Vector{Int}, k::Int, p::Int)
    @show cols
    @assert length(cols) == k
    n, kp = size(data)
    @assert kp == k * p

    @assert size(factors, 1) == n
    n_taxa = div(size(factors, 2), k)
    @assert n_taxa * k == size(factors, 2)

    for i = 1:k
        offset = (i - 1) * p
        check_ind = offset + cols[i]
        fac_inds = i:k:(n_taxa * k)
        for j = 1:n
            if data[j, check_ind] < 0.0
                for k = (offset + 1):(i * p)
                    data[j, k] = -data[j, k]
                    loadings[j, k] = -loadings[j, k]
                end
                for ind in fac_inds
                    factors[j, ind] = -factors[j, ind]
                end
            end
        end
    end
    return data
end


function svd_logs(path::String, new_path::String, k::Int, p::Int;
        cols::Vector{Int} = zeros(Int, k),
        Lid::String = "L",
        fid::String = "factors.",
        rotate_factors::Bool = false,
        transposed::Bool = false)
    d, v, l, f = process_log(path, Lid, k, p, rotation_cols = cols,
                                transposed = transposed,
                                fac_header = fid,
                                rotate_factors = rotate_factors)
    d_labels = ["sv$i" for i = 1:k]
    v_labels = vec(["V$i$j" for j = 1:p, i = 1:k])
    l_labels = vec(["L$i$j" for j = 1:p, i = 1:k])
    cols = get_cols(path)
    f_inds = rotate_factors ? findall(x -> startswith(x, fid), cols) : Int[]
    f_labels = string.(cols[f_inds])


    states = get_log(path, "state", burnin = 0.0)[2]

    labels = ["state"; d_labels; v_labels; l_labels; f_labels]
    data = [states d v l f]


    make_log(new_path, data, labels, includes_states = true)
end

function svd_logs(path::String, new_path::String;
    prec_start = "factorPrecision",
    Lid::String = "L",
    fid::String = "factors.",
    rotate_factors::Bool = false,
    transposed::Bool = false)

    cols = Logs.get_cols(path)
    p = length(findall(x -> startswith(x, prec_start), cols))
    kp = length(findall(x -> startswith(x, Lid), cols))
    k, r = divrem(kp, p)
    if r != 0
        error("Unable to determine the number of factors and traits.")
    end

    return svd_logs(path, new_path, k, p, Lid = Lid, fid = fid,
            rotate_factors = rotate_factors, transposed = transposed)
end



# k = 2
# p = 100
#
# path = joinpath(Directories.desktop, "simFactor_N100_P100_K2_gibbs.log")
# new_path = joinpath(Directories.desktop, "gibbs.log")
#
# svd_logs(path, new_path, k, p)
end
