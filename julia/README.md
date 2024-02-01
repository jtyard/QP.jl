# Julia workflow
The [julia package manager](https://github.com/JuliaLang/juliaup) is probably the easiest way to get julia.  

For generic linux systems on Intel x86 (running natively, in [Windows 11 WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) and on a [Chromebook](https://chromeos.dev/en/linux), you can also install julia with the following commands from a bash shell
```
$ wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.0-linux-x86_64.tar.gz
tar -xvf julia-1.10.0-linux-x86_64.tar.gz
sudo rm /usr/local/bin/julia  # only if previously installed
sudo ln -s ./julia-1.10.0/bin/julia /usr/local/bin/julia
sudo chmod +x /usr/local/bin/julia # may not be needed but can't hurt
``` 

I don't recomment installing julia directly in Windows.  I don't think this is supported by Oscar. 

I keep these packages [`OhMyREPL`](https://kristofferc.github.io/OhMyREPL.jl/latest/), 
[`PkgTemplates`](https://github.com/JuliaCI/PkgTemplates.jl),
[`Revise`](https://github.com/timholy/Revise.jl),
[`Term`](https://fedeclaudi.github.io/Term.jl/stable/),
[`Pkg`](https://github.com/JuliaLang/Pkg.jl) installed in my default environment.  From a fresh install this can be done by pressing `]` to get to the pkg shell then doing 
```
pkg> add OhMyREPL, PkgTemplates, Revise, Term, Pkg
```

My [`startup.jl`](startup.jl) is set up to automatically load an enviroment if julia is started from a directory containing a `Project.toml` or `Manifest.toml`.  It also loads some helpful things for navigating types.  It goes in `.julia/config/startup.jl`.  

To get a local copy of the code to a project you can do e.g.

```
pkg> dev Oscar
```

and this will put it in `.julia/dev/Oscar/`.  You can update the package by running `git pull` from inside `.julia/dev/Oscar/` and you can also navigate that code by adding `.julia/dev` to your vscode workspace.  

