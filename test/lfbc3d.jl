using FBCPoisson
using Test
using LinearAlgebra
using Random
using NonuniformFFTs
using SpecialFunctions: erf

@inline function gaussian_laplace3d_pot(center::NTuple{3, Float64}, target::NTuple{3, Float64}, sigma::Float64)
    r = norm(target .- center)
    return erf(r / (sqrt(2.0) * sigma)) / r / (4π)
end

@testset "NUFFT (m, sigma) selection vs high-accuracy reference" begin
    Random.seed!(2026)
    Ns = (32, 32, 32)
    M = 200
    points = (rand(M) .* 2π, rand(M) .* 2π, rand(M) .* 2π)
    nonuniform_data = randn(ComplexF64, M)
    uniform_data = randn(ComplexF64, Ns...)

    pref = PlanNUFFT(ComplexF64, Ns; m = HalfSupport(15), σ = 2.0)
    set_points!(pref, points)
    type1_ref = Array{ComplexF64}(undef, size(pref))
    type2_ref = Vector{ComplexF64}(undef, M)
    exec_type1!(type1_ref, pref, nonuniform_data)
    exec_type2!(type2_ref, pref, uniform_data)

    prev_m = 0
    for tol in (1e-3, 1e-6, 1e-9)
        sigma_sel, m_sel = FBCPoisson._select_nufft_sigma_m(tol; dim = 3)
        @test m_sel >= prev_m
        prev_m = m_sel

        psel = PlanNUFFT(ComplexF64, Ns; m = HalfSupport(m_sel), σ = sigma_sel)
        set_points!(psel, points)

        type1_sel = Array{ComplexF64}(undef, size(psel))
        type2_sel = Vector{ComplexF64}(undef, M)
        exec_type1!(type1_sel, psel, nonuniform_data)
        exec_type2!(type2_sel, psel, uniform_data)

        relerr1 = norm(type1_sel .- type1_ref) / norm(type1_ref)
        relerr2 = norm(type2_sel .- type2_ref) / norm(type2_ref)
        @test relerr1 <= tol
        @test relerr2 <= tol
    end
end

@inline function gaussian_laplace3d_grad(center::NTuple{3, Float64}, target::NTuple{3, Float64}, sigma::Float64)
    r_vec = target .- center
    r = norm(r_vec)
    a = 1 / (sqrt(2.0) * sigma)
    ddr = ((2a / sqrt(π)) * exp(-(a * r)^2) * r - erf(a * r)) / (r^2)
    return (ddr / (4π * r)) .* collect(r_vec)
end

@testset "lfbc3d vs Gaussian analytic solution" begin
    # Match BoundaryIntegral.jl gaussian_laplace3d_{pot,grad}: a normalized Gaussian density.
    center = (0.03, -0.02, 0.01)
    sigma = 0.08

    n = 48
    xs = collect(range(-0.5, 0.5; length = n))
    h = xs[2] - xs[1]
    w = h^3
    ns = n^3
    sources = Matrix{Float64}(undef, 3, ns)
    charges = Vector{Float64}(undef, ns)
    idx = 1
    for x in xs, y in xs, z in xs
        sources[:, idx] .= (x, y, z)
        r2 = (x - center[1])^2 + (y - center[2])^2 + (z - center[3])^2
        rho = exp(-r2 / (2 * sigma^2)) / ((2π)^(3 / 2) * sigma^3)
        charges[idx] = rho * w
        idx += 1
    end

    targets = [
         0.33 -0.11  0.07 -0.26 0.01
        -0.15  0.20 -0.29  0.31 0.01
         0.10 -0.27  0.24 -0.08 0.01
    ]

    exact_pot = [
        gaussian_laplace3d_pot(center, (targets[1, i], targets[2, i], targets[3, i]), sigma)
        for i in axes(targets, 2)
    ]
    exact_grad = hcat([
        gaussian_laplace3d_grad(center, (targets[1, i], targets[2, i], targets[3, i]), sigma)
        for i in axes(targets, 2)
    ]...)

    N = 64
    pot = lfbc3d(N, sources, charges, targets, 1)
    @test size(pot) == (size(targets, 2),)
    @test maximum(abs.(pot .- exact_pot)) < 6e-3

    pot2, grad2 = lfbc3d(N, sources, charges, targets, 2)
    @test maximum(abs.(pot2 .- exact_pot)) < 6e-3
    @test size(grad2) == size(exact_grad)
    @test maximum(abs.(grad2 .- exact_grad)) < 3e-2

    pre = lfbc3d_precompute(N, sources, charges)
    phase = lfbc3d_prepare_evaluation(pre, targets)
    pot3, grad3 = lfbc3d_evaluate(pre, phase, 2)
    @test maximum(abs.(pot3 .- pot2)) < 1e-12
    @test maximum(abs.(grad3 .- grad2)) < 1e-12
    @test pre.m >= 2
    @test pre.sigma in (1.25, 2.0)

    pre2 = lfbc3d_precompute(N, sources, charges; sigma = 2.0, m = 6)
    @test pre2.sigma == 2.0
    @test pre2.m == 6
end
