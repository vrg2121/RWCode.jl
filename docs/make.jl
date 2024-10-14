push!(LOAD_PATH,"..")
using Documenter, RWCode

makedocs(
    modules = [RWCode],
    sitename="RWCode.jl",
    authors="Conor Walsh and Victoria Garcia",
    pages=["Home" => "index.md"]
)

deploydocs(repo = "github.com/vrg2121/RWCode.jl")