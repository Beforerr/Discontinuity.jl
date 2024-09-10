default:
    just --list

install:
    julia --project -e 'using Pkg; Pkg.develop(["Beforerr", "PlasmaFormulary"]); Pkg.instantiate()'