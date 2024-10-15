push!(LOAD_PATH, "../src/")

using Documenter, RWCode

makedocs(sitename = "RWCode.jl Documentation")

deplydocs(
    repo = "github.com/vrg2121/RWCode.jl.git"
)