struct LinearInterpolationDimension{tType <: AbstractVector} <:
       AbstractInterpolationDimension
    t::tType
    iguesser::Guesser{tType}
    linear_lookup::Bool
    function LinearInterpolationDimension(t::tType; assume_linear_t = 1e-2) where {tType}
        linear_lookup = seems_linear(assume_linear_t, t)
        new{tType}(t, Guesser(t), linear_lookup)
    end
end
