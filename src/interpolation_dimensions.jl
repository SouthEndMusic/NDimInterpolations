struct LinearInterpolationDimension{
    tType <: AbstractVector,
    t_evalType <: AbstractVector,
    idxType <: AbstractVector
} <: AbstractInterpolationDimension
    t::tType
    t_eval::t_evalType
    idx_eval::idxType
    function LinearInterpolationDimension(t, t_eval, idx_eval)
        new{typeof(t), typeof(t_eval), typeof(idx_eval)}(t, t_eval, idx_eval)
    end
end

@adapt_structure LinearInterpolationDimension

function LinearInterpolationDimension(
        t::tType;
        t_eval::t_evalType = similar(t, 0)
) where {tType, t_evalType}
    idx_eval = similar(t_eval, Int)
    itp_dim = LinearInterpolationDimension(
        t, t_eval, idx_eval
    )
    set_eval_idx!(itp_dim)
    itp_dim
end

struct ConstantInterpolationDimension{
    tType <: AbstractVector,
    t_evalType <: AbstractVector,
    idxType <: AbstractVector
} <: AbstractInterpolationDimension
    t::tType
    left::Bool
    t_eval::t_evalType
    idx_eval::idxType
    function ConstantInterpolationDimension(t, left, t_eval, idx_eval)
        new{typeof(t), typeof(t_eval), typeof(idx_eval)}(
            t, left, t_eval, idx_eval
        )
    end
end

@adapt_structure ConstantInterpolationDimension

function ConstantInterpolationDimension(
        t::tType;
        left::Bool = true,
        t_eval::t_evalType = similar(t, 0)
) where {tType, t_evalType}
    idx_eval = similar(t_eval, Int)
    itp_dim = ConstantInterpolationDimension(
        t, left, t_eval, idx_eval
    )
    set_eval_idx!(itp_dim)
    itp_dim
end