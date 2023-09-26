#https://github.com/JuliaLang/Pkg.jl
print("Loading Pkg.jl")
@time import Pkg

# https://github.com/JuliaCI/PkgTemplates.jl 
print("Loading PkgTemplates.jl")
@time import PkgTemplates

# https://github.com/timholy/Revise.jl
print("Loading Revise.jl ")
@time using Revise   # automatically track changes 

# https://kristofferc.github.io/OhMyREPL.jl/latest/
print("Loading OhMyREPL.jl ")
@time using OhMyREPL # syntax highlighting etc in the REPL

# https://fedeclaudi.github.io/Term.jl/stable/
print("Loading Term.jl ")
@time import Term: typestree, termshow, Introspection.inspect 

function create_package(package_name) 
    t = package_template()
    t(package_name)
    println("Created package")
end

function package_template()
    PkgTemplates.Template(; 
        user="yourgithubusername",
        authors="Your Name",
        julia=v"1.9.3",
        dir="~",
        plugins=[
            PkgTemplates.License(; name="GPL-3.0+"),
            # Don't commit Manifest.jl, use SSH and append .jl to the name of the repo
            # See https://invenia.github.io/PkgTemplates.jl/stable/user/#PkgTemplates.Git for more options
            PkgTemplates.Git(; manifest=false, ssh=true, jl=true),
            # Add generated packages to the current environment by `dev`ing them
            # See https://julialang.github.io/Pkg.jl/v1/managing-packages/#Developing-packages-1
            
            #Develop(),
        ]
    )
end

println("\n Use create_package(\"MyPackage\") to create MyPackage.jl in ~/MyPackage/ \n Use inspect(MyType), termshow(MyType), typestree(MyType) to examine types and their dependencies.\n")

if isfile("Project.toml") || isfile("Manifest.toml")
    Pkg.activate(".")
end

