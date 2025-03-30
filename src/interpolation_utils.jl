trivial_range(i::Integer) = i:i

Base.length(itp_dim::AbstractInterpolationDimension) = length(itp_dim.t)

function validate_derivative_orders(derivative_orders::NTuple{N_in, <:Integer}) where {N_in}
    @assert all(≥(0), derivative_orders)
end

function get_ts(interpolation_dimensions::NTuple{
        N_in, AbstractInterpolationDimension}) where {N_in}
    ntuple(i -> interpolation_dimensions[i].t, N_in)
end

function collect_caches(interp::NDimInterpolation{N_out, N_in}) where {N_out, N_in}
    (; u, interpolation_dimensions) = interp
    ts = get_ts(interpolation_dimensions)
    t_evals = ntuple(i -> interpolation_dimensions[i].t_eval, N_in)
    idx_evals = ntuple(i -> interpolation_dimensions[i].idx_eval, N_in)
    u, ts, t_evals, idx_evals
end

function get_output_size(interp::NDimInterpolation{N_out, N_in}) where {N_out, N_in}
    size(interp.u)[(N_in + 1):end]
end

# TODO: Implement a more efficient (GPU compatible) version
function get_idx(tvec::AbstractVector{<:Number}, t_eval::Number)
    idx = 0
    for t in tvec
        if t_eval < t
            break
        else
            idx += 1
        end
    end
    clamp(idx, 1, length(tvec) - 1)
end

function get_idx(
        t::Tuple{Vararg{Number, N_in}}, interpolation_dimensions::NTuple{N_in}) where {N_in}
    idx_eval = ntuple(dim_in -> begin
            itp_dim = interpolation_dimensions[dim_in]
            get_idx(itp_dim.t, t[dim_in])
        end, N_in)
    return idx_eval
end

function get_idx!(
        idx::AbstractVector{Int},
        tvec::AbstractVector{<:Number},
        ts_eval::AbstractVector{<:Number}
)
    backend = get_backend(tvec)
    if !isempty(idx)
        get_idx_kernel(backend)(
            idx,
            tvec,
            ts_eval,
            ndrange = length(ts_eval)
        )
    end
    synchronize(backend)
end

@kernel function get_idx_kernel(
        idx,
        @Const(tvec),
        @Const(t_eval)
)
    i = @index(Global, Linear)
    idx[i] = get_idx(tvec, t_eval[i])
end