struct LinearInterpolationDimension{
    tType <: AbstractVector,
    t_evalType <: AbstractVector,
    idxType <: AbstractVector
} <: AbstractInterpolationDimension
    t::tType
    t_eval::t_evalType
    idx_eval::idxType
    function LinearInterpolationDimension(
            t::tType;
            t_eval::t_evalType = similar(t, 0)
    ) where {tType, t_evalType}
        idx_eval = similar(t_eval, Int)
        get_idx!(idx_eval, t, t_eval)
        new{tType, t_evalType, typeof(idx_eval)}(
            t, t_eval, idx_eval
        )
    end
end
