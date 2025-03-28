using Documenter, NDimInterpolations

makedocs(
    sitename = "NDimInterpolations.jl",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/DataInterpolations/stable/"
    ),
    pages = ["index.md", "usage.md"]
)
