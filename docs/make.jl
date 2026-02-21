using Documenter, TrendDecomposition

makedocs(
    sitename = "TrendDecomposition",
    modules = [TrendDecomposition],
format = Documenter.HTML(
    canonical = "https://sdbrinkmann.github.io/TrendDecomposition.jl/stable/",
    footer = "© 2026 Stefan D. Brinkmann",
    ),
    pages = [
	"Introduction" => "index.md",
	"Get Started" => "man/start.md",
        "Moving Average" => "man/moving.md",
        "Penalized Smoothing" => "man/penalized.md",
        "Exponential Smoothing" => "man/exponential.md",
        "Spectral Analysis" => "man/spectral.md",
	"API" => [
            "Miscellaneous" => "man/misc.md",
            "Index" => "man/api.md",
        ],
    ]	
)


deploydocs(
    repo = "github.com/sdBrinkmann/TrendDecomposition.jl.git",
)
