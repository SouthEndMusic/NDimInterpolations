function _interpolate!(
        out,
        A::NDInterpolation{N_out, N_in, <:LinearInterpolationDimension},
        t::Tuple{Vararg{Number, N_in}},
        idx::NTuple{N_in, <:Integer},
        derivative_orders::NTuple{N_in, <:Integer}
) where {N_out, N_in}
    if iszero(N_out)
        out = zero(out)
    else
        out .= 0
    end
    any(>(1), derivative_orders) && return out

    tᵢ = ntuple(i -> A.interp_dims[i].t[idx[i]], N_in)
    tᵢ₊₁ = ntuple(i -> A.interp_dims[i].t[idx[i] + 1], N_in)

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
        J = (ntuple(i -> idx[i] + I[i], N_in)..., ntuple(_ -> Colon(), N_out)...)
        if iszero(N_out)
            out += c * A.u[J...]
        else
            @. out += c * A.u[J...]
        end
    end
    return out
end

function _interpolate!(
        out,
        A::NDInterpolation{N_out, N_in, <:ConstantInterpolationDimension},
        t::Tuple{Vararg{Number, N_in}},
        idx::NTuple{N_in, <:Integer},
        derivative_orders::NTuple{N_in, <:Integer}
) where {N_out, N_in}
    if iszero(N_out)
        out = zero(out)
    else
        out .= 0
    end
    any(>(0), derivative_orders) && return out
    if iszero(N_out)
        out = A.u[idx...]
    else
        out .= A.u[idx...]
    end
    return out
end