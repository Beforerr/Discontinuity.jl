default:
    just --list

install:
   #!/usr/bin/env julia --project
   using Pkg;
   Pkg.develop(url="https://github.com/Beforerr/PlasmaFormulary.jl");
   Pkg.develop(url="https://github.com/Beforerr/beforerr.jl");
   Pkg.instantiate();