using Documenter
using DarkIntegers


makedocs(
    modules = [DarkIntegers],
    format = :html,
    sitename = "DarkIntegers.jl",
    authors = "Bogdan Opanchuk",
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "API reference" => "api.md",
        "Version history" => "history.md",
    ],
    html_prettyurls = false,
)

deploydocs(
    repo = "github.com/nucypher/DarkIntegers.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
