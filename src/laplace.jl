function truncated_laplace3d(x::T, L::T) where T
    return x > L ? zero(T) : 1 / (4 * π * x)
end

# k is the wavenumber, L is the truncation length
function truncated_laplace3d_hat(k::T, L::T) where T
    if iszero(k)
        return L^2 / 2
    end
    return 2 * (sin(L * k / 2) / k)^2
end
