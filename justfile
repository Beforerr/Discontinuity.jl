default:
    just --list

readme:
    quarto render index.qmd -o README.md -t gfm
    cp README.md docs/src/index.md
