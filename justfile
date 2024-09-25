default:
    just --list

install-dev-deps:
    #!/usr/bin/env -S julia --project
    using Pkg;
    Beforerr = PackageSpec(url="https://github.com/Beforerr/beforerr.jl");
    PlasmaFormulary = PackageSpec(url="https://github.com/Beforerr/PlasmaFormulary.jl");
    Pkg.develop([Beforerr, PlasmaFormulary]);

install: install-dev-deps
    #!/usr/bin/env  julia --project
    using Pkg;
    Pkg.instantiate();
