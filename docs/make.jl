#using Pkg
#Pkg.activate("..")
#push!(LOAD_PATH,"../src/")
using Documenter, RWCode

makedocs(sitename = "RWCode.jl",
    pages = [
        "Home" => "index.md",
        "About" => "about.md",
        "Quick Start" => "quickstart.md",
        "Installation" => "installation.md",
        "Tutorial" => "tutorial.md",
        "Model Configuration" => "configuration.md",
        "Function: run_rwcode()" => "runmodel.md"
    ]
)

deploydocs(
    repo = "github.com/vrg2121/RWCode.jl.git",
    versions = nothing
)