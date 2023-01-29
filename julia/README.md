# Julia workflow



I keep these packages [`OhMyREPL`](https://kristofferc.github.io/OhMyREPL.jl/latest/), 
[`PkgTemplates`](https://github.com/JuliaCI/PkgTemplates.jl),
[`Revise`](https://github.com/timholy/Revise.jl),
[`Term`](https://fedeclaudi.github.io/Term.jl/stable/),
[`Pkg`](https://github.com/JuliaLang/Pkg.jl) installed in my default environment.  From a fresh install this can be done by pressing `]` to get to the pkg shell then doing 
```
pkg> add OhMyREPL, PkgTemplates, Revise, Term, Pkg
```

My [`startup.jl`](startup.jl) is set up to automatically load an enviroment if julia is started from a directory containing a `Project.toml` and `Manifest.toml`.  It also loads some helpful things for navigating types.  It goes in `.julia/config/startup.jl`.  

To get a local copy of the code to a project do e.g. 

```
pkg> dev Oscar
```

and this will put it in `.julia/dev/Oscar/`.  You can update the package by running `git pull` from inside `.julia/dev/Oscar/` and you can also navigate that code e.g. by adding `.julia/dev` to your vscode workspace.  

