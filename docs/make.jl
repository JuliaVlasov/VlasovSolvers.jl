push!(LOAD_PATH,"../src/")

using Documenter
using Literate
using Plots # to not capture precompilation output
using VlasovSolvers

ENV["GKSwstype"] = "100"

examples = [
"Vlasov-Ampere"  => "vlasov-ampere",
"Vlasov-Poisson" => "vlasov-poisson",
"Bump On Tail"   => "bump_on_tail",
"Rotation 2D"    => "rotation2d_bsl",
"Vlasov-HMF"     => "vlasov-hmf"
]

# generate examples

for example in examples

    @show EXAMPLE     = joinpath(@__DIR__, "examples", string(example[2],".jl"))
    @show DOC_OUTPUT  = joinpath(@__DIR__, "src")
    @show NB_OUTPUT   = joinpath(@__DIR__, "src", "notebooks")
   
    Literate.markdown(EXAMPLE, DOC_OUTPUT)
    Literate.notebook(EXAMPLE, NB_OUTPUT, execute=false)

end

makedocs(modules=[VlasovSolvers],
         sitename = "VlasovSovers.jl",
         authors="Pierre Navaro",
         format=Documenter.HTML(;
         prettyurls=get(ENV, "CI", "false") == "true",
         canonical="https://juliavlasov.github.io/VlasovSolvers.jl",
         assets=String[],
         ),
         doctest = false,
         pages = ["Home"     => "index.md",
		  "Contents" => "contents.md",
                  "Examples" => [string(example[2],".md") for example in examples]])

deploydocs(;
    branch = "gh-pages",
    devbranch = "master",
    repo   = "github.com/JuliaVlasov/VlasovSolvers.jl"
)
