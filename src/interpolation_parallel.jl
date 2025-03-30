function eval_unstructured(
        interp::NDimInterpolation{N_out, N_in};
        kwargs...
) where {N_out, N_in}
    n_points = length(first(interp.interpolation_dimensions).t_eval)
    out = similar(interp.u, (n_points, get_output_size(interp)...))
    eval_unstructured!(out, interp; kwargs...)
end

function eval_unstructured!(
        out::AbstractArray,
        interp::NDimInterpolation{N_out, N_in, I},
        derivative_orders::NTuple{N_in, <:Integer} = ntuple(_ -> 0, N_in)
) where {N_out, N_in, I}
    backend = get_backend(out)
    u, ts, t_evals, idx_evals = collect_caches(interp)
    @assert all(t -> length(t) == size(out, 1), t_evals)
    @assert size(out)[2:end] == get_output_size(interp)
    eval_kernel(backend)(
        out,
        u, ts, t_evals, idx_evals,
        derivative_orders,
        :unstructured,
        N_out, N_in, I,
        ndrange = size(out, 1)
    )
    synchronize(backend)
    return out
end

function eval_grid(interp::NDimInterpolation{N_out, N_in}; kwargs...) where {N_out, N_in}
    grid_size = ntuple(i -> interp.interpolation_dimensions[i].t_eval, N_in)
    out = similar(interp.u, (grid_size..., get_output_size(interp)...))
    eval_grid!(out, interp; kwargs...)
end

function eval_grid!(
        out::AbstractArray,
        interp::NDimInterpolation{N_out, N_in, I};
        derivative_orders::NTuple{N_in, <:Integer} = ntuple(_ -> 0, N_in)
) where {N_out, N_in, I}
    backend = get_backend(out)
    u, ts, t_evals, idx_evals = collect_caches(interp)
    @assert all(i -> size(out, i) == length(t_evals[i]), N_in)
    @assert size(out)[(N_in + 1):end] == get_output_size(interp)
    eval_kernel(backend)(
        out,
        u, ts, t_evals, idx_evals,
        derivative_orders,
        :grid,
        N_out, N_in, I,
        ndrange = size(out)[1:N_in]
    )
    synchronize(backend)
    return out
end

@kernel function eval_kernel(
        out,
        @Const(u),
        @Const(ts),
        @Const(t_evals),
        @Const(idx_evals),
        derivative_orders,
        mode,
        N_out,
        N_in,
        I
)
    k = @index(Global, NTuple)

    if mode == :unstructured
        t_eval = ntuple(i -> t_evals[i][only(k)], N_in)
        idx_eval = ntuple(i -> idx_evals[i][only(k)], N_in)
    else # mode == :grid
        t_eval = ntuple(i -> t_evals[i][k[i]], N_in)
        idx_eval = ntuple(i -> idx_evals[i][k[i]], N_in)
    end

    if iszero(N_out)
        out[k...] = _interpolate(u, ts, t_eval, idx_eval, derivative_orders, I)
    else
        _interpolate!(
            view(out, k..., ntuple(_ -> Colon(), N_out)...),
            u, ts, t_eval, idx_eval, derivative_orders, I)
    end
end