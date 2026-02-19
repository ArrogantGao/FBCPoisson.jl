using FBCPoisson
using Documenter

DocMeta.setdocmeta!(FBCPoisson, :DocTestSetup, :(using FBCPoisson); recursive=true)

makedocs(;
    modules=[FBCPoisson],
    authors="Xuanzhao Gao <xgao@flatironinstitute.org> and contributors",
    sitename="FBCPoisson.jl",
    format=Documenter.HTML(;
        canonical="https://ArrogantGao.github.io/FBCPoisson.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ArrogantGao/FBCPoisson.jl",
    devbranch="main",
)
