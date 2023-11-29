push!(LOAD_PATH,"../src/")

using Documenter
using Literate
using Plots 
using VlasovSolvers

ENV["GKSwstype"] = "100"

makedocs(modules=[VlasovSolvers],
         sitename = "VlasovSolvers.jl",
         authors="Julia Vlasov",
         format=Documenter.HTML(;
         prettyurls=get(ENV, "CI", "false") == "true",
         canonical="https://juliavlasov.github.io/VlasovSolvers.jl",
         assets=String[],
         ),
         doctest = false,
         pages = [ "Home"          => "index.md",
                   "Quickstart"    => "problem.md",
                   "Vlasov-Ampere" => "vlasov-ampere.md",
                   "Bump On Tail"  => "bump_on_tail.md",
                   "Rotation 2D"   => "rotation2d.md",
                   "Vlasov-HMF"    => "vlasov-hmf.md",
                   "Two-stream instability" => "tsi.md",
                   "API"           => "api.md",
                   "Contents"      => "contents.md"])

deploydocs(;
    branch = "gh-pages",
    devbranch = "master",
    repo   = "github.com/JuliaVlasov/VlasovSolvers.jl"
)
