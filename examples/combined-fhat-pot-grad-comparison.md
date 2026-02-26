# Combined Comparison: Fourier Mode Error and Pot/Grad Errors

- FINUFFT backend
- Bandwidths: (0.03, 0.08, 0.18, 0.30, 0.5, 1.0)
- Tolerances: (1e-3, 1e-6, 1e-9)
- Setup: N=128, L=1.8, source span (default), nsrc=40 uniform, targets=10x10x10 in [-0.5,0.5]^3
- Note: If source span is too small for the target box, pot/grad columns are reported as N/A.

Columns:
1) relerr_fhat_type1_vs_analytic
2) relerr_pot_analytic_fhat, relerr_grad_analytic_fhat
3) relerr_pot_type1_fhat, relerr_grad_type1_fhat

| bw | tol | relerr_fhat_type1_vs_analytic | relerr_pot_analytic_fhat | relerr_grad_analytic_fhat | relerr_pot_type1_fhat | relerr_grad_type1_fhat |
|---:|---:|---:|---:|---:|---:|---:|
| 0.03 | 1e-03 | 3.573e-04 | 5.767e-04 | 8.396e-04 | 5.461e-04 | 8.089e-04 |
| 0.03 | 1e-06 | 2.255e-07 | 7.678e-05 | 7.993e-04 | 7.678e-05 | 7.993e-04 |
| 0.03 | 1e-09 | 2.846e-11 | 7.678e-05 | 7.994e-04 | 7.678e-05 | 7.994e-04 |
| 0.08 | 1e-03 | 1.063e-03 | 2.585e-04 | 2.423e-04 | 3.487e-04 | 3.156e-04 |
| 0.08 | 1e-06 | 3.944e-08 | 9.469e-09 | 2.083e-08 | 1.192e-08 | 2.461e-08 |
| 0.08 | 1e-09 | 9.291e-10 | 5.683e-10 | 7.157e-10 | 5.987e-10 | 6.813e-10 |
| 0.18 | 1e-03 | 1.063e-03 | 1.784e-04 | 2.640e-04 | 3.465e-04 | 3.118e-04 |
| 0.18 | 1e-06 | 3.944e-08 | 1.792e-08 | 3.198e-08 | 1.813e-08 | 3.773e-08 |
| 0.18 | 1e-09 | 9.291e-10 | 5.133e-10 | 6.196e-10 | 5.925e-10 | 6.244e-10 |
| 0.30 | 1e-03 | 1.063e-03 | 1.964e-04 | 2.317e-04 | 1.919e-04 | 2.574e-04 |
| 0.30 | 1e-06 | 3.944e-08 | 1.275e-08 | 1.555e-08 | 1.337e-08 | 1.996e-08 |
| 0.30 | 1e-09 | 9.291e-10 | 5.155e-10 | 6.088e-10 | 5.909e-10 | 5.881e-10 |
| 0.50 | 1e-03 | 1.063e-03 | 1.435e-04 | 3.202e-04 | 2.136e-04 | 5.008e-04 |
| 0.50 | 1e-06 | 3.944e-08 | 5.016e-09 | 2.215e-08 | 6.685e-09 | 2.275e-08 |
| 0.50 | 1e-09 | 9.291e-10 | 4.246e-10 | 6.603e-10 | 3.927e-10 | 5.664e-10 |
| 1.00 | 1e-03 | 1.063e-03 | 1.948e-04 | 4.832e-04 | 2.064e-04 | 6.981e-04 |
| 1.00 | 1e-06 | 3.944e-08 | 6.957e-09 | 3.569e-08 | 8.101e-09 | 3.971e-08 |
| 1.00 | 1e-09 | 9.291e-10 | 4.681e-10 | 1.147e-09 | 4.536e-10 | 1.085e-09 |
