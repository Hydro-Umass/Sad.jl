push!(LOAD_PATH, "../src")

import Pkg
Pkg.activate("..")

using Documenter, Sad

makedocs(sitename="Sad.jl")
