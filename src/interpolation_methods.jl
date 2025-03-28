# Out of place KA friendly scalar valued interpolation
function _interpolate(
        u::AbstractArray{Tu, N_in},
        ts::NTuple{N_in, AbstractVector{Tt}},
        t::Tuple{Vararg{Number, N_in}},
        idxs::NTuple{N_in, <:Integer},
        derivative_orders::NTuple{N_in, <:Integer},
        ::Type{<:LinearInterpolationDimension} # For dispatching, can probably be done better
) where {Tu, N_in, Tt}
    out = zero(promote_type(Tu, Tt, map(typeof, t)...))
    any(>(1), derivative_orders) && return out

    tᵢ = ntuple(i -> ts[i][idxs[i]], N_in)
    tᵢ₊₁ = ntuple(i -> ts[i][idxs[i] + 1], N_in)

    # Size of the (hyper)volume `t` is in
    t_vol = one(eltype(tᵢ))
    for (t₁, t₂) in zip(tᵢ, tᵢ₊₁)
        t_vol *= t₂ - t₁
    end

    # Loop over corners of the interval product `t` is in
    for I in Iterators.product(ntuple(i -> (false, true), N_in)...)
        c = eltype(out)(inv(t_vol))
        for (t_, right_point, d, t₁, t₂) in zip(t, I, derivative_orders, tᵢ, tᵢ₊₁)
            c *= if right_point
                iszero(d) ? t_ - t₁ : one(t_)
            else
                iszero(d) ? t₂ - t_ : -one(t_)
            end
        end
        J = ntuple(i -> idxs[i] + I[i], N_in)
        out += c * u[J...]
    end
    return out
end

# In-place KA friendly vector valued interpolation
function _interpolate!(
        out,
        u::AbstractArray{T, N},
        ts::NTuple{N_in, <:AbstractArray},
        t::Tuple{Vararg{Number, N_in}},
        idxs::NTuple{N_in, <:Integer},
        derivative_orders::NTuple{N_in, <:Integer},
        ::Type{<:LinearInterpolationDimension} # For dispatching, can probably be done better
) where {T, N, N_in}
    out .= 0
    any(>(1), derivative_orders) && return out

    tᵢ = ntuple(i -> ts[i][idxs[i]], N_in)
    tᵢ₊₁ = ntuple(i -> ts[i][idxs[i] + 1], N_in)

    # Size of the (hyper)volume `t` is in
    t_vol = one(eltype(tᵢ))
    for (t₁, t₂) in zip(tᵢ, tᵢ₊₁)
        t_vol *= t₂ - t₁
    end

    # Loop over corners of the interval product `t` is in
    for I in Iterators.product(ntuple(i -> (false, true), N_in)...)
        c = eltype(out)(inv(t_vol))
        for (t_, right_point, d, t₁, t₂) in zip(t, I, derivative_orders, tᵢ, tᵢ₊₁)
            c *= if right_point
                iszero(d) ? t_ - t₁ : one(t_)
            else
                iszero(d) ? t₂ - t_ : -one(t_)
            end
        end
        for J in CartesianIndices(
            ntuple(i -> i ≤ N_in ? trivial_range(idxs[i] + I[i]) : 1:size(u, i), N),
        )
            out[Tuple(J)[(N_in + 1):N]...] += c * u[J]
        end
    end
    return out
end
