trivial_range(i::Integer) = i:i

Base.length(itp_dim::AbstractInterpolationDimension) = length(itp_dim.t)

function validate_derivative_orders(derivative_orders::NTuple{N_in, <:Integer}) where {N_in}
    @assert all(≥(0), derivative_orders)
end

function get_ts(interp_dims::NTuple{
        N_in, AbstractInterpolationDimension}) where {N_in}
    ntuple(i -> interp_dims[i].t, N_in)
end

function get_output_size(interp::NDInterpolation{N_out, N_in}) where {N_out, N_in}
    size(interp.u)[(N_in + 1):end]
end

function make_out(
        interp::NDInterpolation{0, N_in},
        t::NTuple{N_in, >:Number}
) where {N_in}
    zero(promote_type(eltype(interp.u), map(typeof, t)...))
end

function make_out(
        interp::NDInterpolation{N_out, N_in},
        t::NTuple{N_in, >:Number}
) where {N_out, N_in}
    similar(
        interp.u, promote_type(eltype(interp.u), map(eltype, t)...), get_output_size(interp))
end

# TODO: Implement a more efficient (GPU compatible) version
function get_idx(
        interp_dim::AbstractInterpolationDimension,
        t_eval::Number
)
    (; t) = interp_dim
    left = (interp_dim isa ConstantInterpolationDimension) ? interp_dim.left : true
    idx = clamp(searchsortedlast(t, t_eval), 1, length(t) - 1)
    !left && (idx += 1)
    idx
end

function get_idx(
        interp_dims::NTuple{N_in},
        t::Tuple{Vararg{Number, N_in}};
) where {N_in}
    ntuple(dim_in -> get_idx(interp_dims[dim_in], t[dim_in]), N_in)
end

function set_eval_idx!(
        interp_dim::AbstractInterpolationDimension,
)
    backend = get_backend(interp_dim.t)
    if !isempty(interp_dim.t_eval)
        set_idx_kernel(backend)(
            interp_dim,
            ndrange = length(interp_dim.t_eval)
        )
    end
    synchronize(backend)
end

@kernel function set_idx_kernel(
        interp_dim
)
    i = @index(Global, Linear)
    interp_dim.idx_eval[i] = get_idx(interp_dim, interp_dim.t_eval[i])
end