push!(LOAD_PATH, "../src")

import Pkg
Pkg.add("Documenter")
Pkg.activate("..")

using Documenter, Sad

makedocs(sitename="Sad.jl",
         checkdocs = :none,
         format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"))
