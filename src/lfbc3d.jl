function _as_3xn(points::AbstractMatrix{T}) where {T <: Real}
    if size(points, 1) == 3
        return Matrix{T}(points)
    elseif size(points, 2) == 3
        return permutedims(points, (2, 1))
    else
        throw(ArgumentError("points must have shape (3, N) or (N, 3)"))
    end
end

struct LFBC3DPrecomputation{T, A3}
    Fk::A3
    k1::Vector{T}
    center::NTuple{3, T}
    scale::T
    Δk::T
    nufft_tol::T
    quad_weight::T
end

struct LFBC3DEvaluationPhase{T, V, A3}
    nt::Int
    pot_complex::V
    grad_complex::V
    Fg::A3
    targets_unit::Matrix{T}
    trg_points::NTuple{3, Vector{T}}
end

"""
    lfbc3d_precompute(N, sources, targets, charges, nufft_tol)

Precompute source-dependent data for 3D free-boundary Poisson evaluation:
- selects normalization box from combined source/target extents
- computes type-1 NUFFT coefficients
- applies truncated Green's function in Fourier space.
`L` and `margin` are fixed internally.
"""
function lfbc3d_precompute(
    N::Integer,
    sources::AbstractMatrix{<:Real},
    targets::AbstractMatrix{<:Real},
    charges::AbstractVector{<:Real},
    nufft_tol::Real,
)
    L = 1.8
    margin = 1.0

    N > 0 || throw(ArgumentError("N must be positive"))

    T = promote_type(Float64, eltype(sources), eltype(targets), eltype(charges), typeof(L), typeof(margin), typeof(nufft_tol))
    src = _as_3xn(Matrix{T}(sources))
    q = Vector{T}(charges)
    trg = _as_3xn(Matrix{T}(targets))
    size(src, 2) == length(q) || throw(DimensionMismatch("number of charges must match number of sources"))

    mins = min.(minimum(src; dims = 2), minimum(trg; dims = 2))
    maxs = max.(maximum(src; dims = 2), maximum(trg; dims = 2))

    center_vec = vec((mins .+ maxs) ./ 2)
    center = (center_vec[1], center_vec[2], center_vec[3])

    span = maximum(maxs .- mins)
    scale = iszero(span) ? one(T) : span / T(margin)

    src_unit = (src .- reshape(collect(center), 3, 1)) ./ scale

    # Δk = T(π / 2)
    Δk = T(2π / (1 + sqrt(3)))
    src_points = (Δk .* vec(src_unit[1, :]), Δk .* vec(src_unit[2, :]), Δk .* vec(src_unit[3, :]))
    Nint = Int(N)
    Fk = nufft3d1(src_points[1], src_points[2], src_points[3], complex.(q), -1, T(nufft_tol), Nint, Nint, Nint)

    # FINUFFT uses centered mode ordering: -N/2, ..., N/2-1 (for even N).
    k1 = Vector{T}(undef, Nint)
    for i in 1:Nint
        m = i - 1 - (Nint ÷ 2)
        k1[i] = Δk * m
    end

    Lt = T(L)
    for i in 1:Nint, j in 1:Nint, k in 1:Nint
        kval = sqrt(k1[i]^2 + k1[j]^2 + k1[k]^2)
        Fk[i, j, k] *= truncated_laplace3d_hat(kval, Lt)
    end

    quad_weight = (Δk / (2π))^3
    return LFBC3DPrecomputation(Fk, k1, center, scale, Δk, T(nufft_tol), quad_weight)
end

"""
    lfbc3d_prepare_evaluation(pre, targets)

Prepare target-dependent evaluation phase (target points, buffers and phase state).
"""
function lfbc3d_prepare_evaluation(pre::LFBC3DPrecomputation{T}, targets::AbstractMatrix{<:Real}) where {T}
    trg = _as_3xn(Matrix{T}(targets))
    center = reshape(collect(pre.center), 3, 1)
    targets_unit = (trg .- center) ./ pre.scale

    if any(abs.(targets_unit) .>= T(0.5))
        throw(ArgumentError("normalized targets hit periodic boundary; increase precompute span or margin"))
    end

    trg_points = NTuple{3, Vector{T}}((
        pre.Δk .* vec(targets_unit[1, :]),
        pre.Δk .* vec(targets_unit[2, :]),
        pre.Δk .* vec(targets_unit[3, :]),
    ))

    nt = size(trg, 2)
    pot_complex = Vector{Complex{T}}(undef, nt)
    grad_complex = Vector{Complex{T}}(undef, nt)
    Fg = similar(pre.Fk)
    return LFBC3DEvaluationPhase(nt, pot_complex, grad_complex, Fg, targets_unit, trg_points)
end

"""
    lfbc3d_evaluate(pre, phase, pgt)

Evaluate potential (`pgt = 1`) or potential+gradient (`pgt = 2`) using precomputed data.
"""
function lfbc3d_evaluate(pre::LFBC3DPrecomputation{T}, phase::LFBC3DEvaluationPhase{T}, pgt::Integer) where {T}
    pgt in (1, 2) || throw(ArgumentError("pgt must be 1 (potential) or 2 (potential + gradient)"))

    phase.pot_complex .= nufft3d2(
        phase.trg_points[1],
        phase.trg_points[2],
        phase.trg_points[3],
        1,
        pre.nufft_tol,
        pre.Fk,
    )
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
        phase.grad_complex .= nufft3d2(
            phase.trg_points[1],
            phase.trg_points[2],
            phase.trg_points[3],
            1,
            pre.nufft_tol,
            phase.Fg,
        )
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
    nufft_tol::Real,
    pgt::Integer;
)
    pre = lfbc3d_precompute(N, sources, targets, charges, nufft_tol)
    phase = lfbc3d_prepare_evaluation(pre, targets)
    return lfbc3d_evaluate(pre, phase, pgt)
end
