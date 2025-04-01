function eval_unstructured(
        interp::NDInterpolation{N_out, N_in};
        kwargs...
) where {N_out, N_in}
    n_points = length(first(interp.interp_dims).t_eval)
    out = similar(interp.u, (n_points, get_output_size(interp)...))
    eval_unstructured!(out, interp; kwargs...)
end

function eval_unstructured!(
        out::AbstractArray,
        interp::NDInterpolation{N_out, N_in, I},
        derivative_orders::NTuple{N_in, <:Integer} = ntuple(_ -> 0, N_in)
) where {N_out, N_in, I}
    validate_derivative_orders(derivative_orders)
    backend = get_backend(out)
    @assert all(i -> length(interp.interp_dims[i].t_eval) == size(out, 1), N_in)
    @assert size(out)[2:end] == get_output_size(interp)
    eval_kernel(backend)(
        out,
        interp,
        derivative_orders,
        false,
        ndrange = size(out, 1)
    )
    synchronize(backend)
    return out
end

function eval_grid(interp::NDInterpolation{N_out, N_in}; kwargs...) where {N_out, N_in}
    grid_size = ntuple(i -> interp.interp_dims[i].t_eval, N_in)
    out = similar(interp.u, (grid_size..., get_output_size(interp)...))
    eval_grid!(out, interp; kwargs...)
end

function eval_grid!(
        out::AbstractArray,
        interp::NDInterpolation{N_out, N_in, I};
        derivative_orders::NTuple{N_in, <:Integer} = ntuple(_ -> 0, N_in)
) where {N_out, N_in, I}
    validate_derivative_orders(derivative_orders)
    backend = get_backend(out)
    @assert all(i -> size(out, i) == length(interp.interp_dims[i].t_eval), N_in)
    @assert size(out)[(N_in + 1):end] == get_output_size(interp)
    eval_kernel(backend)(
        out,
        interp,
        derivative_orders,
        true,
        ndrange = size(out)[1:N_in]
    )
    synchronize(backend)
    return out
end

@kernel function eval_kernel(
        out,
        @Const(A),
        derivative_orders,
        eval_grid
)
    N_in = length(A.interp_dims)
    N_out = ndims(A.u) - N_in

    k = @index(Global, NTuple)

    if eval_grid
        t_eval = ntuple(i -> A.interp_dims[i].t_eval[k[i]], N_in)
        idx_eval = ntuple(i -> A.interp_dims[i].idx_eval[k[i]], N_in)
    else
        t_eval = ntuple(i -> A.interp_dims[i].t_eval[only(k)], N_in)
        idx_eval = ntuple(i -> A.interp_dims[i].idx_eval[only(k)], N_in)
    end

    if iszero(N_out)
        out[k...] = _interpolate!(
            make_out(A, t_eval), A, t_eval, idx_eval, derivative_orders)
    else
        _interpolate!(
            view(out, k..., ntuple(_ -> Colon(), N_out)...),
            A, t_eval, idx_eval, derivative_orders)
    end
end