using Documenter, NDInterpolations

makedocs(
    sitename = "NDInterpolations.jl",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/DataInterpolations/stable/"
    ),
    pages = [
        "index.md",
        "Usage" => "usage.md",
        "Interpolation Types" => "interpolation_types.md"
    ]
)
