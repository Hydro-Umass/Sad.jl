import Pkg
Pkg.add("Documenter")

push!(LOAD_PATH, "../src")

using Documenter, Sad

makedocs(
    modules = [Sad],
    sitename="Sad.jl",
    checkdocs = :none,
         format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"))
