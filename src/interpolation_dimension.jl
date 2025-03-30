struct LinearInterpolationDimension{
    tType <: AbstractVector,
    t_evalType <: AbstractVector,
    idxType <: AbstractVector
} <: AbstractInterpolationDimension
    t::tType
    iguesser::Guesser{tType}
    linear_lookup::Bool
    t_eval::t_evalType
    idx_eval::idxType
    function LinearInterpolationDimension(
            t::tType;
            assume_linear_t = 1e-2,
            t_eval::t_evalType = similar(t, 0)
    ) where {tType, t_evalType}
        linear_lookup = seems_linear(assume_linear_t, t)
        iguesser = Guesser(t)
        idx_eval = get_idx.(Ref(t), t_eval, Ref(iguesser))
        new{tType, t_evalType, typeof(idx_eval)}(
            t, iguesser, linear_lookup, t_eval, idx_eval)
    end
end
