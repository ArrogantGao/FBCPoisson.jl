using FBCPoisson
using Test
using FFTW

@testset "FBCPoisson.jl" begin
    include("laplace.jl")
    include("lfbc3d.jl")
end
