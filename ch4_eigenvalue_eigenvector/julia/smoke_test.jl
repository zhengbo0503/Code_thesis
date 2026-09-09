# smoke_test.jl - Chapter 4 Julia smoke test
#
# Exercises the vendored JacobiEigen package: the plain Jacobi eigensolver
# and both mixed-precision variants (mp2, mp3), checking orthogonality,
# reconstruction, eigenvalue ordering, and that mp3 attains better relative
# accuracy on tiny eigenvalues than the plain Float64 solver.
#
# Run from ch4_eigenvalue_eigenvector/julia with:
#     julia --project=JacobiEigen smoke_test.jl

using LinearAlgebra, JacobiEigen, Quadmath, Printf

npass = 0; nfail = 0
function chk(name, cond, detail="")
    global npass, nfail
    cond ? (npass += 1) : (nfail += 1)
    @printf("  [%s] %-44s %s\n", cond ? "PASS" : "FAIL", name, detail)
end

println("=== Chapter 4 (Julia): JacobiEigen ===")
println("  Julia ", VERSION)

# SPD test matrix with a controlled spread of eigenvalues
n = 40
using Random; Random.seed!(0)
Q, _ = qr(randn(n, n)); Q = Matrix(Q)
d = 10.0 .^ range(0, -10, length = n)
A = Q * Diagonal(d) * Q'; A = (A + A') / 2
dref = sort(Float64.(eigvals(Float128.(A))), rev = true)   # high-precision reference

# --- off --------------------------------------------------------------
# NOTE: JacobiEigen.off computes sqrt(norm(A)^2 - sum|A_ii|^2), which cancels
# catastrophically. On a diagonal matrix the true value is 0 but it returns
# ~norm(A)*sqrt(eps) ~ 1.8e-8, so off cannot resolve anything below that
# floor. MATLAB's off.m zeroes the diagonal and takes the Frobenius norm
# instead, and is exact. Asserted here at the achievable floor, not at 0.
offdiag = JacobiEigen.off(Matrix(Diagonal(d)))
chk("off(diagonal) at cancellation floor", offdiag < 1e-7,
    @sprintf("%.3e (exact answer is 0; sqrt-eps floor)", offdiag))
chk("off(A) > 0", JacobiEigen.off(A) > 0, @sprintf("%.3e", JacobiEigen.off(A)))

# --- plain Jacobi -------------------------------------------------------
Λ, V, P = jacobi_eigen(A)
chk("jacobi_eigen V orthogonal", norm(V'V - I) < 1e-12, @sprintf("%.3e", norm(V'V - I)))
chk("jacobi_eigen reconstruction", norm(A - V*Diagonal(Λ)*V') / norm(A) < 1e-12,
    @sprintf("%.3e", norm(A - V*Diagonal(Λ)*V') / norm(A)))
chk("jacobi_eigen eigenvalues sorted", issorted(Λ) || issorted(Λ, rev = true))
chk("jacobi_eigen params (rot, sweep)", length(P) == 2 && P[2] < 30, string(P))

# --- mixed precision ----------------------------------------------------
Λ2, V2, P2, t2 = mp2_jacobi_eigen(copy(A), Float32)
chk("mp2 V orthogonal", norm(V2'V2 - I) < 1e-12, @sprintf("%.3e", norm(V2'V2 - I)))
chk("mp2 reconstruction", norm(A - V2*Diagonal(Λ2)*V2') / norm(A) < 1e-12,
    @sprintf("%.3e", norm(A - V2*Diagonal(Λ2)*V2') / norm(A)))
chk("mp2 timing vector length 4", length(t2) == 4, @sprintf("%.3fs total", sum(t2)))

Λ3, V3, P3, t3 = mp3_jacobi_eigen(copy(A), Float32, Float128)
chk("mp3 V orthogonal", norm(V3'V3 - I) < 1e-12, @sprintf("%.3e", norm(V3'V3 - I)))
chk("mp3 reconstruction", norm(A - V3*Diagonal(Λ3)*V3') / norm(A) < 1e-12,
    @sprintf("%.3e", norm(A - V3*Diagonal(Λ3)*V3') / norm(A)))
chk("mp3 timing vector length 4", length(t3) == 4, @sprintf("%.3fs total", sum(t3)))

# --- the accuracy claim -------------------------------------------------
relerr(x) = maximum(abs.(sort(Float64.(x), rev = true) .- dref) ./ abs.(dref))
e_lapack = relerr(eigvals(A))
e_jac    = relerr(Λ)
e_mp3    = relerr(Λ3)
@printf("  rel err: LAPACK=%.3e  Jacobi=%.3e  mp3=%.3e\n", e_lapack, e_jac, e_mp3)
chk("mp3 more accurate than LAPACK eigvals", e_mp3 < e_lapack,
    @sprintf("mp3=%.3e lapack=%.3e", e_mp3, e_lapack))

@printf("  ---- Chapter 4 Julia: %d/%d passed, %d failed ----\n", npass, npass + nfail, nfail)
exit(nfail == 0 ? 0 : 1)
