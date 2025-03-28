module NDimInterpolations
import FindFirstFunctions: searchsortedfirstcorrelated, searchsortedlastcorrelated,
                           Guesser

abstract type AbstractInterpolationDimension end

include("interpolation_utils.jl")
include("interpolation_dimension.jl")
include("interpolation_methods.jl")

struct NDimInterpolation{
    N_out, N_in, I <: AbstractInterpolationDimension, uType <: AbstractArray}
    interpolation_dimensions::NTuple{N_in, I}
    u::uType
    function NDimInterpolation(
            interpolation_dimensions::NTuple{N_in, I},
            u::AbstractArray{T, N}
    ) where {N_in, I, T, N}
        N_out = N - N_in # Compile time inferrable?
        @assert N_out ≥ 0
        @assert ntuple(i -> length(interpolation_dimensions[i]), N_in) == size(u)[1:N_in]
        new{N_out, N_in, I, typeof(u)}(interpolation_dimensions, u)
    end
end

# Multiple `t` arguments to tuple (can these 2 be done in 1?)
function (interp::NDimInterpolation)(t_args::Vararg{<:Number}; kwargs...)
    interp(t_args; kwargs...)
end

function (interp::NDimInterpolation)(
        out::AbstractArray, t_args::Vararg{<:Number}; kwargs...)
    interp(out, t_args; kwargs...)
end

# Out of place single input evaluation
function (interp::NDimInterpolation{
        N_out, N_in, I})(
        t::Tuple{Vararg{Number, N_in}};
        derivative_orders::NTuple{N_in, <:Integer} = ntuple(_ -> 0, N_in)
) where {N_out, N_in, I, tType <: Number}
    validate_derivative_orders(derivative_orders)
    (; interpolation_dimensions, u) = interp
    idxs, ts = get_inputs(t, interpolation_dimensions)
    if iszero(N_out)
        # Scalar output
        _interpolate(u, ts, t, idxs, derivative_orders, I)
    else
        # Vector output
        out = similar(
            u, promote_type(eltype(u), map(eltype, t)...), size(u)[(N_in + 1):end])
        _interpolate!(out, u, ts, t, idxs, derivative_orders, I)
    end
end

# In place single input evaluation
function (interp::NDimInterpolation{N_out, N_in, I})(
        out::AbstractArray,
        t::Tuple{Vararg{Number, N_in}};
        derivative_orders::NTuple{N_in, <:Integer} = ntuple(_ -> 0, N_in)
) where {N_out, N_in, I}
    validate_derivative_orders(derivative_orders)
    (; interpolation_dimensions, u) = interp
    idxs, ts = get_inputs(t, interpolation_dimensions)
    @assert size(out) == size(u)[(N_in + 1):end]
    _interpolate!(out, u, ts, t, idxs, derivative_orders, I)
end

export NDimInterpolation, LinearInterpolationDimension

end # module NDimInterpolations
