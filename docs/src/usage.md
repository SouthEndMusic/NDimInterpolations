## Single point evaluation

### Setup

```@example tutorial
using NDimInterpolations
using Random
Random.seed!(1)

interpolation_dimensions = (
    LinearInterpolationDimension(cumsum(0.5 .+ rand(5))),
    LinearInterpolationDimension(cumsum(0.5 .+ rand(10)))
)

t_eval_1 = range(
    first(interpolation_dimensions[1].t),
    last(interpolation_dimensions[1].t),
    length = 100
)
t_eval_2 = range(
    first(interpolation_dimensions[2].t),
    last(interpolation_dimensions[2].t),
    length = 100
)
```

### Scalar output

```@example tutorial
using Plots

u = rand(5, 10)
itp = NDimInterpolation(interpolation_dimensions, u)
out = itp.(t_eval_1, t_eval_2')
heatmap(out)
```

### Vector output, out of place

```@example tutorial
u = reshape(u, 5, 10, 1)
itp = NDimInterpolation(interpolation_dimensions, u)
out = itp.(t_eval_1, t_eval_2')
heatmap(map(only, out))
```

### vector output, in place

```@example tutorial
out = zeros(100, 100)
for I in CartesianIndices(out)
    i, j = Tuple(I)
    itp(view(out, i, j:j), t_eval_1[i], t_eval_2[j])
end
heatmap(out)
```

## Single point derivative evaluation

### Partial derivative w.r.t. first input (ForwardDiff)

```@example tutorial
using ForwardDiff
u = reshape(u, 5, 10)
itp = NDimInterpolation(interpolation_dimensions, u)

```

### Partial derivative w.r.t. first input (analytic)

```@example tutorial

```