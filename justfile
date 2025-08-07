default:
    just --list

readme:
    quarto render index.qmd -o README.md -t gfm
    cp README.md docs/src/index.md


install:
    julia -e '
    using Pkg
    Pkg.develop([
        PackageSpec(url="https://github.com/JuliaPlasma/PlasmaFormulary.jl"),
        PackageSpec(url="https://github.com/Beforerr/Beforerr.jl"),
        PackageSpec(path=pwd())
    ])
    Pkg.instantiate()
    '