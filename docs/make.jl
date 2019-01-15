using Documenter
using DarkIntegers


makedocs(
    modules = [DarkIntegers],
    format = Documenter.HTML(prettyurls=false),
    sitename = "DarkIntegers.jl",
    authors = "Bogdan Opanchuk",
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "API reference" => "api.md",
        "Version history" => "history.md",
    ],
)

deploydocs(
    repo = "github.com/nucypher/DarkIntegers.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
