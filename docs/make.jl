using Documenter
using ReedMuller

DocMeta.setdocmeta!(ReedMuller, :DocTestSetup, :(using ReedMuller); recursive = true)

makedocs(;
    modules = [ReedMuller],
    sitename = "ReedMuller.jl",
    authors = "Anatoly Rodionov",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://ayarodionov.github.io/ReedMuller.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/ayarodionov/ReedMuller.jl",
    devbranch = "main",
)
