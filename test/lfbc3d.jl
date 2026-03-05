using FBCPoisson
using Test
using LinearAlgebra
using Random
using FINUFFT
using SpecialFunctions: erf

@inline function gaussian_laplace3d_pot(center::NTuple{3, Float64}, target::NTuple{3, Float64}, sigma::Float64)
    r = norm(target .- center)
    return erf(r / (sqrt(2.0) * sigma)) / r / (4π)
end

@inline function gaussian_laplace3d_grad(center::NTuple{3, Float64}, target::NTuple{3, Float64}, sigma::Float64)
    r_vec = target .- center
    r = norm(r_vec)
    a = 1 / (sqrt(2.0) * sigma)
    ddr = ((2a / sqrt(π)) * exp(-(a * r)^2) * r - erf(a * r)) / (r^2)
    return (ddr / (4π * r)) .* collect(r_vec)
end

function gauss_legendre_nodes_weights(n::Int, a::Float64, b::Float64)
    @assert n >= 2
    β = [k / sqrt(4k^2 - 1) for k in 1:(n - 1)]
    T = SymTridiagonal(zeros(n), β)
    F = eigen(T)
    x = collect(F.values)
    w = collect(2 .* (F.vectors[1, :]).^2)
    x = ((b - a) / 2) .* x .+ (a + b) / 2
    w = ((b - a) / 2) .* w
    return x, w
end

function make_source_quadrature(kind::Symbol, n::Int, region::Float64, center::NTuple{3, Float64}, sigma::Float64)
    if kind == :uniform
        xs = collect(range(-region, region; length = n))
        h = xs[2] - xs[1]
        wx = fill(h, n)
        ys, wy = xs, wx
        zs, wz = xs, wx
    elseif kind == :gl
        xs, wx = gauss_legendre_nodes_weights(n, -region, region)
        ys, wy = gauss_legendre_nodes_weights(n, -region, region)
        zs, wz = gauss_legendre_nodes_weights(n, -region, region)
    else
        error("unknown quadrature kind: $kind")
    end

    ns = n^3
    sources = Matrix{Float64}(undef, 3, ns)
    charges = Vector{Float64}(undef, ns)
    idx = 1
    for ix in eachindex(xs), iy in eachindex(ys), iz in eachindex(zs)
        x = xs[ix]
        y = ys[iy]
        z = zs[iz]
        sources[:, idx] .= (x, y, z)
        r2 = (x - center[1])^2 + (y - center[2])^2 + (z - center[3])^2
        rho = exp(-r2 / (2 * sigma^2)) / ((2π)^(3 / 2) * sigma^3)
        charges[idx] = rho * wx[ix] * wy[iy] * wz[iz]
        idx += 1
    end
    return sources, charges
end

function make_targets_10x10x10(a::Float64, b::Float64)
    xs = collect(range(a, b; length = 10))
    nt = 10^3
    targets = Matrix{Float64}(undef, 3, nt)
    idx = 1
    for x in xs, y in xs, z in xs
        targets[:, idx] .= (x, y, z)
        idx += 1
    end
    return targets
end

@inline gaussian_fhat_abs(k::Float64, bw::Float64) = exp(-0.5 * (bw * k)^2)

function select_N_from_gaussian_fhat(
    bw::Float64,
    tol::Float64;
    Δk::Float64 = π / 2,
    Nmin::Int = 64,
    safety::Float64 = 2.0,
)
    @assert bw > 0
    @assert tol > 0
    @assert safety >= 1
    kreq = sqrt(2 * log(1 / tol)) / bw
    N = max(Nmin, ceil(Int, safety * kreq / Δk))
    return iseven(N) ? N : N + 1
end

@testset "NUFFT tol vs high-accuracy reference" begin
    Random.seed!(2026)
    Ns = (32, 32, 32)
    M = 200
    points = (rand(M) .* 2π, rand(M) .* 2π, rand(M) .* 2π)
    nonuniform_data = randn(ComplexF64, M)
    uniform_data = randn(ComplexF64, Ns...)

    eps_ref = 1e-15
    type1_ref = nufft3d1(points[1], points[2], points[3], nonuniform_data, 1, eps_ref, Ns...)
    type2_ref = nufft3d2(points[1], points[2], points[3], 1, eps_ref, uniform_data)

    for tol in (1e-3, 1e-6, 1e-9)
        type1_sel = nufft3d1(points[1], points[2], points[3], nonuniform_data, 1, tol, Ns...)
        type2_sel = nufft3d2(points[1], points[2], points[3], 1, tol, uniform_data)

        relerr1 = norm(type1_sel .- type1_ref) / norm(type1_ref)
        relerr2 = norm(type2_sel .- type2_ref) / norm(type2_ref)
        @test relerr1 <= 5tol
        @test relerr2 <= 5tol
    end
end

@testset "lfbc3d pot/grad: bandwidth x tol x quadrature" begin
    center = (0.1, - 0.2, 0.3)
    bandwidths = (0.2, 0.5, 1.0, 2.0)
    tols = (1e-3, 1e-6, 1e-9)
    quadratures = (:uniform,)

    targets = make_targets_10x10x10(-0.6, 0.4)
    for bw in bandwidths
        # Choose domain half-width R so exp(-R^2 / (2*bw^2)) = 1e-12.
        region = bw * sqrt(2 * log(1e12))
        nsrc = 40
        sources, charges = make_source_quadrature(:uniform, nsrc, region, center, bw)

        # Analytic Gaussian reference.
        exact_pot = [
            gaussian_laplace3d_pot(center, (targets[1, i], targets[2, i], targets[3, i]), bw)
            for i in axes(targets, 2)
        ]
        exact_grad = hcat([
            gaussian_laplace3d_grad(center, (targets[1, i], targets[2, i], targets[3, i]), bw)
            for i in axes(targets, 2)
        ]...)

        prev_N = 0
        for tol in tols

            @info "Testing bw=$bw, tol=$tol."

            # N = select_N_from_gaussian_fhat(bw, tol)
            # @test gaussian_fhat_abs(N * (π / 2), bw) <= tol
            # @test N >= prev_N
            N = 128
            prev_N = N

            pre = lfbc3d_precompute(N, sources, targets, charges, tol)

            phase = lfbc3d_prepare_evaluation(pre, targets)
            pot, grad = lfbc3d_evaluate(pre, phase, 2)

            relerr_pot = norm(pot .- exact_pot) / norm(exact_pot)
            relerr_grad = norm(grad .- exact_grad) / norm(exact_grad)
            @test relerr_pot <= tol * 5
            @test relerr_grad <= tol * 5
        end

    end
end
