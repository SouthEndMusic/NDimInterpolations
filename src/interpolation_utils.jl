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

##
### NOTE: All functions below are copied from DataInterpolations.jl but have to be replaced by KA friendly versions
##

seems_linear(assume_linear_t::Bool, _) = assume_linear_t
seems_linear(assume_linear_t::Number, t) = looks_linear(t; threshold = assume_linear_t)

function looks_linear(t; threshold = 1e-2)
    length(t) <= 2 && return true
    t_0, t_f = first(t), last(t)
    t_span = t_f - t_0
    tspan_over_N = t_span * length(t)^(-1)
    norm_var = sum(
        (t_i - t_0 - i * tspan_over_N)^2 for (i, t_i) in enumerate(t)
    ) / (length(t) * t_span^2)
    norm_var < threshold^2
end

function get_idx(tvec::AbstractVector{<:Number}, t_eval::Number,
        iguess::Union{<:Integer, Guesser}; lb = 1,
        ub_shift = -1, idx_shift = 0, side = :last)
    ub = length(tvec) + ub_shift
    return if side == :last
        clamp(searchsortedlastcorrelated(tvec, t_eval, iguess) + idx_shift, lb, ub)
    elseif side == :first
        clamp(searchsortedfirstcorrelated(tvec, t_eval, iguess) + idx_shift, lb, ub)
    else
        error("side must be :first or :last")
    end
end

function get_idx(
        t::Tuple{Vararg{Number, N_in}}, interpolation_dimensions::NTuple{N_in}) where {N_in}
    idx_eval = ntuple(dim_in -> begin
            itp_dim = interpolation_dimensions[dim_in]
            get_idx(itp_dim.t, t[dim_in], itp_dim.iguesser)
        end, N_in)
    return idx_eval
end