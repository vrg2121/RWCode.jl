push!(LOAD_PATH, "../src/")

using Documenter, RWCode

makedocs(sitename = "RWCode.jl Documentation",
    format = Documenter.HTML(),
    source = "src/"
)

deploydocs(
    repo = "github.com/vrg2121/RWCode.jl.git",
    target = "build",
    branch = "gh-pages"
)