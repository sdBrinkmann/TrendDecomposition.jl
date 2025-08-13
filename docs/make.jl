using Documenter, TrendDecomposition

makedocs(
    sitename = "TrendDecomposition",
    modules = [TrendDecomposition],
    pages = [
	"Introduction" => "index.md",
	"Get Started" => "man/start.md",
        "Moving Average" => "man/moving.md",
        "Penalized Smoothing" => "man/penalized.md",
        "Exponential Smoothing" => "man/exponential.md",
	"API" => "man/api.md",
    ]	
)


deploydocs(
    repo = "github.com/sdBrinkmann/TrendDecomposition.jl.git",
)
