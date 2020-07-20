import Pkg

packages = ["Blink", "Interact", "CSV", "DataFrames", "RCall", "Statistics"]

Pkg.add(packages)
Pkg.add(Pkg.PackageSpec(url="https://github.com/gabehassler/BeastUtils.jl.git"))
