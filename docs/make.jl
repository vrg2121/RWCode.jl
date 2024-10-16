push!(LOAD_PATH, "../src/")

using Documenter, RWCode

makedocs(sitename = "RWCode.jl Documentation"
)

deploydocs(
    repo = "github.com/vrg2121/RWCode.jl.git",
    branch = "gh-pages"
)