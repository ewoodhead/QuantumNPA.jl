using Documenter
using QuantumNPA

makedocs(sitename="QuantumNPA.jl",
         authors = "Erik Woodhead",
         pages = ["Home" => "index.md",
                  "Examples" => "examples.md",
                  "NPA" => "npa.md"])

deploydocs(repo = "github.com/ewoodhead/QuantumNPA.jl.git")
