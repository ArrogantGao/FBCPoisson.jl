using FBCPoisson
using Test
using FFTW

@testset "truncated Laplace kernel FFT consistency (3D)" begin
    L = 1.8
    N = 64
    box = 4L
    dx = box / N

    x = ((0:(N - 1)) .- N ÷ 2) .* dx

    g = Array{Float64}(undef, N, N, N)
    for ix in 1:N, iy in 1:N, iz in 1:N
        r = sqrt(x[ix]^2 + x[iy]^2 + x[iz]^2)
        # The kernel is singular at r=0 but integrable in 3D.
        # Setting this single grid point to 0 introduces only high-order quadrature error.
        g[ix, iy, iz] = iszero(r) ? 0.0 : truncated_laplace3d(r, L)
    end

    Gnum = dx^3 .* fftshift(fft(ifftshift(g)))
    k = ((0:(N - 1)) .- N ÷ 2) .* (2π / box)

    # Compare a few low-frequency modes to limit quadrature/aliasing error.
    modes = (
        (0, 0, 0),
        (1, 0, 0),
        (2, 0, 0),
        (1, 1, 0),
        (2, 1, 0),
        (1, 1, 1),
    )

    for (mx, my, mz) in modes
        ix = mx + N ÷ 2 + 1
        iy = my + N ÷ 2 + 1
        iz = mz + N ÷ 2 + 1

        kabs = sqrt(k[ix]^2 + k[iy]^2 + k[iz]^2)
        expected = truncated_laplace3d_hat(kabs, L)
        observed = real(Gnum[ix, iy, iz])
        @test observed ≈ expected rtol = 3e-2 atol = 2e-3
    end
end
