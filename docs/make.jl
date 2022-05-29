using Documenter, QuadraticOutputSystems

push!(LOAD_PATH, "../src/")
makedocs(
 sitename="QuadraticOutputSystems",
 pages = [
          "Home" => "index.md",
         ]
)
deploydocs(
  repo = "github.com/Algopaul/QuadraticOutputSystems.git",
  versions = nothing
)
