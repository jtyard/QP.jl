
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
    @eval begin
        using PkgTemplates
        Template(; 
            user="jtyard",
            authors="Jon Yard",
            julia=v"1.8",
            dir="~",
            plugins=[
                License(; name="MIT"),
                # Commit Manifest.jl, use SSH and append .jl to the name of the repo
                # See https://invenia.github.io/PkgTemplates.jl/stable/user/#PkgTemplates.Git for more options
                Git(; manifest=true, ssh=true, jl=true),
                # Add generated packages to the current environment by `dev`ing them
                # See https://julialang.github.io/Pkg.jl/v1/managing-packages/#Developing-packages-1
                #Develop(),
            ]
        )(package_name)
    end
end


println("\n Use create_package()(\"MyPackage\") to create MyPackage.jl in ~/MyPackage/ \n Use inspect(MyType), termshow(MyType), typestree(MyType) to examine types and their dependencies.\n")


if isfile("Project.toml") && isfile("Manifest.toml")
    using Pkg
    Pkg.activate(".")
end

