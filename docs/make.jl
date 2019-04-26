using Documenter
using MeshArrays

makedocs(
    sitename = "MeshArrays",
    format = :html,
    modules = [MeshArrays]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
