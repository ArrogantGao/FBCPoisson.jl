function _as_3xn(points::AbstractMatrix{T}) where {T <: Real}
    if size(points, 1) == 3
        return Matrix{T}(points)
    elseif size(points, 2) == 3
        return permutedims(points, (2, 1))
    else
        throw(ArgumentError("points must have shape (3, N) or (N, 3)"))
    end
end

struct LFBC3DPrecomputation{T, P, A3}
    plan::P
    Fk::A3
    k1::Vector{T}
    center::NTuple{3, T}
    scale::T
    Δk::T
    quad_weight::T
    sigma::T
    m::Int
end

struct LFBC3DEvaluationPhase{T, V, A3}
    nt::Int
    pot_complex::V
    grad_complex::V
    Fg::A3
    targets_unit::Matrix{T}
end

function _select_nufft_sigma_m(
    tol::T;
    dim::Int = 3,
    sigma::Union{Nothing, T} = nothing,
    m::Union{Nothing, Int} = nothing,
) where {T <: Real}
    tol_eff = max(tol, Base.eps(T))
    tolfac = T(0.18 * 1.4^(dim - 1))
    nsoff = one(T)

    # Based on FINUFFT discussion #798 (v2.5.0 tuning):
    # m ~= ceil(log(tolfac/tol) / (pi*sqrt(1 - 1/sigma)) + nsoff)
    m_from_sigma(σ::T) = max(2, ceil(Int, log(tolfac / tol_eff) / (π * sqrt(1 - 1 / σ)) + nsoff))

    if !isnothing(m) && !isnothing(sigma)
        return T(sigma), Int(m)
    elseif !isnothing(m)
        # If only m is given, keep robust default sigma.
        return T(2.0), Int(m)
    elseif !isnothing(sigma)
        return T(sigma), m_from_sigma(T(sigma))
    else
        # For stricter tolerances, prefer sigma = 2.0 for robustness.
        if tol_eff <= T(1e-5)
            σ = T(2.0)
            return σ, m_from_sigma(σ)
        end
        σ_candidates = (T(1.25), T(2.0))
        best = (σ_candidates[1], m_from_sigma(σ_candidates[1]))
        best_cost = best[1]^3 * best[2]
        for σ in σ_candidates[2:end]
            mw = m_from_sigma(σ)
            cost = σ^3 * mw
            if cost < best_cost
                best = (σ, mw)
                best_cost = cost
            end
        end
        return best
    end
end

"""
    lfbc3d_precompute(N, sources, charges; L = 1.8, margin = 0.9, span = nothing)

Precompute source-dependent data for 3D free-boundary Poisson evaluation:
- creates the NUFFT plan
- computes type-1 NUFFT coefficients
- applies truncated Green's function in Fourier space.
"""
function lfbc3d_precompute(
    N::Integer,
    sources::AbstractMatrix{<:Real},
    charges::AbstractVector{<:Real};
    L::Real = 1.8,
    margin::Real = 0.9,
    span::Union{Nothing, Real} = nothing,
    nufft_tol::Real = 1e-6,
    sigma::Union{Nothing, Real} = nothing,
    m::Union{Nothing, Int} = nothing,
)
    N > 0 || throw(ArgumentError("N must be positive"))
    margin > 0 || throw(ArgumentError("margin must be positive"))

    T = promote_type(Float64, eltype(sources), eltype(charges), typeof(L), typeof(margin))
    src = _as_3xn(Matrix{T}(sources))
    q = Vector{T}(charges)
    size(src, 2) == length(q) || throw(DimensionMismatch("number of charges must match number of sources"))

    mins = minimum(src; dims = 2)
    maxs = maximum(src; dims = 2)
    center_vec = vec((mins .+ maxs) ./ 2)
    center = (center_vec[1], center_vec[2], center_vec[3])

    source_span = maximum(maxs .- mins)
    raw_span = isnothing(span) ? source_span : T(span)
    scale = iszero(raw_span) ? one(T) : raw_span / T(margin)

    src_unit = (src .- reshape(collect(center), 3, 1)) ./ scale

    sigma_sel, m_sel = _select_nufft_sigma_m(T(nufft_tol); dim = 3, sigma = isnothing(sigma) ? nothing : T(sigma), m = m)

    Ns = (Int(N), Int(N), Int(N))
    plan = PlanNUFFT(Complex{T}, Ns; m = HalfSupport(m_sel), σ = sigma_sel)

    Δk = T(π / 2)
    src_points = (Δk .* vec(src_unit[1, :]), Δk .* vec(src_unit[2, :]), Δk .* vec(src_unit[3, :]))
    set_points!(plan, src_points)

    Fk = Array{Complex{T}}(undef, size(plan))
    exec_type1!(Fk, plan, complex.(q))

    k1 = Vector{T}(undef, Int(N))
    nhalf = Int(N) ÷ 2
    for i in 1:Int(N)
        m = i <= nhalf ? (i - 1) : (i - 1 - Int(N))
        k1[i] = Δk * m
    end

    Lt = T(L)
    for i in 1:Int(N), j in 1:Int(N), k in 1:Int(N)
        kval = sqrt(k1[i]^2 + k1[j]^2 + k1[k]^2)
        Fk[i, j, k] *= truncated_laplace3d_hat(kval, Lt)
    end

    quad_weight = (Δk / (2π))^3
    return LFBC3DPrecomputation(plan, Fk, k1, center, scale, Δk, quad_weight, sigma_sel, m_sel)
end

"""
    lfbc3d_prepare_evaluation(pre, targets)

Prepare target-dependent evaluation phase (target points, buffers and phase state).
"""
function lfbc3d_prepare_evaluation(
    pre::LFBC3DPrecomputation{T},
    targets::AbstractMatrix{<:Real},
) where {T}
    trg = _as_3xn(Matrix{T}(targets))
    center = reshape(collect(pre.center), 3, 1)
    targets_unit = (trg .- center) ./ pre.scale

    if any(abs.(targets_unit) .>= T(0.5))
        throw(ArgumentError("normalized targets hit periodic boundary; increase precompute span or margin"))
    end

    trg_points = (
        pre.Δk .* vec(targets_unit[1, :]),
        pre.Δk .* vec(targets_unit[2, :]),
        pre.Δk .* vec(targets_unit[3, :]),
    )
    set_points!(pre.plan, trg_points)

    nt = size(trg, 2)
    pot_complex = Vector{Complex{T}}(undef, nt)
    grad_complex = Vector{Complex{T}}(undef, nt)
    Fg = similar(pre.Fk)
    return LFBC3DEvaluationPhase(nt, pot_complex, grad_complex, Fg, targets_unit)
end

"""
    lfbc3d_evaluate(pre, phase, pgt)

Evaluate potential (`pgt = 1`) or potential+gradient (`pgt = 2`) using precomputed data.
"""
function lfbc3d_evaluate(
    pre::LFBC3DPrecomputation{T},
    phase::LFBC3DEvaluationPhase{T},
    pgt::Integer,
) where {T}
    pgt in (1, 2) || throw(ArgumentError("pgt must be 1 (potential) or 2 (potential + gradient)"))

    exec_type2!(phase.pot_complex, pre.plan, pre.Fk)
    pot = pre.quad_weight .* real.(phase.pot_complex) ./ pre.scale

    if pgt == 1
        return pot
    end

    N = length(pre.k1)
    grad = Matrix{T}(undef, 3, phase.nt)
    for d in 1:3
        for i in 1:N, j in 1:N, k in 1:N
            kd = d == 1 ? pre.k1[i] : d == 2 ? pre.k1[j] : pre.k1[k]
            phase.Fg[i, j, k] = (im * kd) * pre.Fk[i, j, k]
        end
        exec_type2!(phase.grad_complex, pre.plan, phase.Fg)
        grad[d, :] .= pre.quad_weight .* real.(phase.grad_complex) ./ (pre.scale^2)
    end

    return pot, grad
end

# solve the potential of a charge distribution under free boundary conditions in 3D
#
# pgt = 1: potential at targets
# pgt = 2: potential and gradient at targets
function lfbc3d(
    N::Integer,
    sources::AbstractMatrix{<:Real},
    charges::AbstractVector{<:Real},
    targets::AbstractMatrix{<:Real},
    pgt::Integer;
    kwargs...,
)
    pre = lfbc3d_precompute(N, sources, charges; kwargs...)
    phase = lfbc3d_prepare_evaluation(pre, targets)
    return lfbc3d_evaluate(pre, phase, pgt)
end
