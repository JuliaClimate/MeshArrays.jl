using Documenter
using MeshArrays

makedocs(
    sitename = "MeshArrays",
    format   = Documenter.HTML(),
    modules  = [MeshArrays],
    pages = [
    "Home" => "index.md",
    "Main" => "main.md",
    "Detail" => "detail.md",
    "API" => "API.md",
            ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/gaelforget/MeshArrays.jl.git",
)
