# Julia workflow



I keep these packages [`OhMyREPL`](https://kristofferc.github.io/OhMyREPL.jl/latest/), 
[PkgTemplates](https://github.com/JuliaCI/PkgTemplates.jl),
[Revise](https://github.com/timholy/Revise.jl),
[Term](https://fedeclaudi.github.io/Term.jl/stable/),
[Pkg](https://github.com/JuliaLang/Pkg.jl) installed in my default environment.  

My [startup.jl](startup file) is set up to automatically load an enviroment if julia is started from a directory containing a `Project.toml` and `Manifest.toml`.  It also loads some helpful things for navigating types.