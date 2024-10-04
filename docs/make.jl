push!(LOAD_PATH,"../src/")
using RWCode
using Documenter
makedocs(
         sitename = "RWCode.jl",
         modules  = [RWCode],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/vrg2121/RWCode.jl",
)