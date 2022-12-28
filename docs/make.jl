import Pkg
Pkg.add("Documenter")

push!(LOAD_PATH, "../src")

using Documenter, Sad

makedocs(
    modules = [Sad],
    sitename="Sad.jl",
    build = "build",
    clean = true,
    doctest = true,
    repo = "https://github.com/Hydro-Umass/Sad.jl",
    checkdocs = :none,
    pages = [
        "Home" => "index.md",
        "Install" => "installation.md",
        "Algorithm" => "algorithm.md",
        "Use cases" => "use_cases.md",
        "API" => "api.md",
    ],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"))

deploydocs(
    repo = "github.com/Hydro-Umass/Sad.jl",
    target = "build",
    versions = nothing,
    push_preview = true,
)
