using FBCPoisson
using FINUFFT
using LinearAlgebra
using SpecialFunctions: erf
using Printf

@inline function gaussian_laplace3d_pot(center::NTuple{3, Float64}, target::NTuple{3, Float64}, sigma::Float64)
    r = norm(target .- center)
    return erf(r / (sqrt(2.0) * sigma)) / r / (4 * pi)
end

@inline function gaussian_laplace3d_grad(center::NTuple{3, Float64}, target::NTuple{3, Float64}, sigma::Float64)
    r_vec = target .- center
    r = norm(r_vec)
    a = 1 / (sqrt(2.0) * sigma)
    ddr = ((2 * a / sqrt(pi)) * exp(-(a * r)^2) * r - erf(a * r)) / (r^2)
    return (ddr / (4 * pi * r)) .* collect(r_vec)
end

function make_source_quadrature_uniform(n::Int, region::Float64, center::NTuple{3, Float64}, sigma::Float64)
    xs = collect(range(-region, region; length = n))
    h = xs[2] - xs[1]
    w = fill(h, n)
    ns = n^3
    sources = Matrix{Float64}(undef, 3, ns)
    charges = Vector{Float64}(undef, ns)
    idx = 1
    for ix in eachindex(xs), iy in eachindex(xs), iz in eachindex(xs)
        x = xs[ix]
        y = xs[iy]
        z = xs[iz]
        sources[:, idx] .= (x, y, z)
        r2 = (x - center[1])^2 + (y - center[2])^2 + (z - center[3])^2
        rho = exp(-r2 / (2 * sigma^2)) / ((2 * pi)^(3 / 2) * sigma^3)
        charges[idx] = rho * w[ix] * w[iy] * w[iz]
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

analytic_fhat_scaled(kx::Float64, ky::Float64, kz::Float64, bw::Float64, scale::Float64) = exp(-0.5 * (bw / scale)^2 * (kx^2 + ky^2 + kz^2))

center = (0.0, 0.0, 0.0)
bandwidths = (0.03, 0.08, 0.18, 0.30, 0.5, 1.0)
tols = (1e-3, 1e-6, 1e-9)
nsrc = 40
N = 128
L = 1.8
targets = make_targets_10x10x10(-0.5, 0.5)

rows = String[]
push!(rows, "| bw | tol | relerr_fhat_type1_vs_analytic | relerr_pot_analytic_fhat | relerr_grad_analytic_fhat | relerr_pot_type1_fhat | relerr_grad_type1_fhat |")
push!(rows, "|---:|---:|---:|---:|---:|---:|---:|")

for bw in bandwidths
    region = bw * sqrt(2 * log(1e12))
    sources, charges = make_source_quadrature_uniform(nsrc, region, center, bw)

    exact_pot = [gaussian_laplace3d_pot(center, (targets[1, i], targets[2, i], targets[3, i]), bw) for i in axes(targets, 2)]
    exact_grad = hcat([gaussian_laplace3d_grad(center, (targets[1, i], targets[2, i], targets[3, i]), bw) for i in axes(targets, 2)]...)

    for tol in tols
        pre = lfbc3d_precompute(N, sources, targets, charges, tol)

        # (1) type-1 NUFFT Fourier mode error vs analytic f_hat
        src_unit = (sources .- reshape(collect(pre.center), 3, 1)) ./ pre.scale
        points = (pre.Δk .* vec(src_unit[1, :]), pre.Δk .* vec(src_unit[2, :]), pre.Δk .* vec(src_unit[3, :]))
        Fk_type1 = nufft3d1(points[1], points[2], points[3], complex.(charges), -1, tol, N, N, N)

        Fk_analytic = Array{ComplexF64}(undef, N, N, N)
        for i in 1:N, j in 1:N, k in 1:N
            Fk_analytic[i, j, k] = analytic_fhat_scaled(pre.k1[i], pre.k1[j], pre.k1[k], bw, pre.scale)
        end
        relerr_fhat = norm(Fk_type1 .- Fk_analytic) / norm(Fk_analytic)

        # (2) pot/grad error using analytic Fourier coeffs
        Fk_analytic_green = similar(Fk_analytic)
        for i in 1:N, j in 1:N, k in 1:N
            kval = sqrt(pre.k1[i]^2 + pre.k1[j]^2 + pre.k1[k]^2)
            Fk_analytic_green[i, j, k] = Fk_analytic[i, j, k] * truncated_laplace3d_hat(kval, L)
        end
        pre_analytic = FBCPoisson.LFBC3DPrecomputation(
            Fk_analytic_green,
            pre.k1,
            pre.center,
            pre.scale,
            pre.Δk,
            pre.nufft_tol,
            pre.quad_weight,
        )
        # (2)+(3) pot/grad comparisons can fail if source-span is too small for the target box.
        relerr_pot_analytic = NaN
        relerr_grad_analytic = NaN
        relerr_pot_type1 = NaN
        relerr_grad_type1 = NaN
        eval_ok = true
        try
            phase_analytic = lfbc3d_prepare_evaluation(pre_analytic, targets)
            pot_analytic, grad_analytic = lfbc3d_evaluate(pre_analytic, phase_analytic, 2)
            relerr_pot_analytic = norm(pot_analytic .- exact_pot) / norm(exact_pot)
            relerr_grad_analytic = norm(grad_analytic .- exact_grad) / norm(exact_grad)

            phase_type1 = lfbc3d_prepare_evaluation(pre, targets)
            pot_type1, grad_type1 = lfbc3d_evaluate(pre, phase_type1, 2)
            relerr_pot_type1 = norm(pot_type1 .- exact_pot) / norm(exact_pot)
            relerr_grad_type1 = norm(grad_type1 .- exact_grad) / norm(exact_grad)
        catch e
            if e isa ArgumentError && occursin("normalized targets hit periodic boundary", sprint(showerror, e))
                eval_ok = false
            else
                rethrow(e)
            end
        end

        if eval_ok
            push!(rows, @sprintf(
                "| %.2f | %.0e | %.3e | %.3e | %.3e | %.3e | %.3e |",
                bw,
                tol,
                relerr_fhat,
                relerr_pot_analytic,
                relerr_grad_analytic,
                relerr_pot_type1,
                relerr_grad_type1,
            ))
        else
            push!(rows, @sprintf(
                "| %.2f | %.0e | %.3e | %s | %s | %s | %s |",
                bw,
                tol,
                relerr_fhat,
                "N/A (target out of box)",
                "N/A (target out of box)",
                "N/A (target out of box)",
                "N/A (target out of box)",
            ))
        end
    end
end

report = joinpath(@__DIR__, "2026-02-26-combined-fhat-pot-grad-comparison.md")
open(report, "w") do io
    println(io, "# Combined Comparison: Fourier Mode Error and Pot/Grad Errors")
    println(io)
    println(io, "- Date: 2026-02-26")
    println(io, "- FINUFFT backend")
    println(io, "- Bandwidths: (0.03, 0.08, 0.18, 0.30)")
    println(io, "- Tolerances: (1e-3, 1e-6, 1e-9)")
    println(io, "- Setup: N=128, L=1.8, source span (default), nsrc=40 uniform, targets=10x10x10 in [-0.5,0.5]^3")
    println(io, "- Note: If source span is too small for the target box, pot/grad columns are reported as N/A.")
    println(io)
    println(io, "Columns:")
    println(io, "1) relerr_fhat_type1_vs_analytic")
    println(io, "2) relerr_pot_analytic_fhat, relerr_grad_analytic_fhat")
    println(io, "3) relerr_pot_type1_fhat, relerr_grad_type1_fhat")
    println(io)
    foreach(r -> println(io, r), rows)
end

println(report)
