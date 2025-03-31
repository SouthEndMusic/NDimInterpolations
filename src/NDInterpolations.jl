module NDInterpolations
using KernelAbstractions # Keep as dependency or make extension?

abstract type AbstractInterpolationDimension end

struct NDInterpolation{
    N_out, N_in, I <: AbstractInterpolationDimension, uType <: AbstractArray}
    interpolation_dimensions::NTuple{N_in, I}
    u::uType
    function NDInterpolation(
            interpolation_dimensions::NTuple{N_in, I},
            u::AbstractArray{T, N}
    ) where {N_in, I, T, N}
        N_out = N - N_in # Compile time inferrable?
        @assert N_out ≥ 0
        @assert ntuple(i -> length(interpolation_dimensions[i]), N_in) == size(u)[1:N_in]
        new{N_out, N_in, I, typeof(u)}(interpolation_dimensions, u)
    end
end

include("interpolation_utils.jl")
include("interpolation_dimension.jl")
include("interpolation_methods.jl")
include("interpolation_parallel.jl")

# Multiple `t` arguments to tuple (can these 2 be done in 1?)
function (interp::NDInterpolation)(t_args::Vararg{<:Number}; kwargs...)
    interp(t_args; kwargs...)
end

function (interp::NDInterpolation)(
        out::AbstractArray, t_args::Vararg{<:Number}; kwargs...)
    interp(out, t_args; kwargs...)
end

# Out of place single input evaluation
function (interp::NDInterpolation{
        N_out, N_in, I})(
        t::Tuple{Vararg{Number, N_in}};
        derivative_orders::NTuple{N_in, <:Integer} = ntuple(_ -> 0, N_in)
) where {N_out, N_in, I}
    validate_derivative_orders(derivative_orders)
    (; interpolation_dimensions, u) = interp
    idx_eval = get_idx(t, interpolation_dimensions)
    ts = get_ts(interpolation_dimensions)
    if iszero(N_out)
        # Scalar output
        _interpolate(u, ts, t, idx_eval, derivative_orders, I)
    else
        # Vector output
        out = similar(
            u, promote_type(eltype(u), map(eltype, t)...), get_output_size(interp))
        _interpolate!(out, u, ts, t, idx_eval, derivative_orders, I)
    end
end

# In place single input evaluation
function (interp::NDInterpolation{N_out, N_in, I})(
        out::AbstractArray,
        t::Tuple{Vararg{Number, N_in}};
        derivative_orders::NTuple{N_in, <:Integer} = ntuple(_ -> 0, N_in)
) where {N_out, N_in, I}
    validate_derivative_orders(derivative_orders)
    (; interpolation_dimensions, u) = interp
    idx_eval = get_idx(t, interpolation_dimensions)
    ts = get_ts(interpolation_dimensions)
    @assert size(out) == size(u)[(N_in + 1):end]
    _interpolate!(out, u, ts, t, idx_eval, derivative_orders, I)
end

export NDInterpolation, LinearInterpolationDimension, eval_unstructured,
       eval_unstructured!, eval_grid, eval_grid!

end # module NDInterpolations
