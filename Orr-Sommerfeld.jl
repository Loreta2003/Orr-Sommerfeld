using LinearAlgebra
using Plots

function fdiff4C(N)
    h = 1 / (N + 1)
    D = diagm(-(N-1) => ones(2), -2 => ones(N-1), -1 => -4*ones(N), 0 => 6*ones(N+1), 1 => -4*ones(N), 2 => ones(N-1), N-1 => ones(2))
    D[end,1] = -4
    D[1,end] = -4
    D = D / h^4
end
function fdiff2C(N)
    h = 1 / (N + 1)
    D = diagm( -1 => ones(N), 0 => -2*ones(N+1), 1 => ones(N))
    D[end,1] = 1
    D[1,end] = 1
    D = D 
end
function orr_sommerfeld_solver(N, Re, alpha, c, U, U_double_prime)
    D2 = fdiff2C(N)
    D4 = fdiff4C(N)

    # Orr-Sommerfeld operator
    A = D4 - 2 * alpha * alpha * D2 + alpha * alpha * alpha * alpha * Diagonal(ones(N+1)) - alpha * Re * 1im * (Diagonal(U) - c*Diagonal(ones(N+1))) * (D2 - alpha * alpha * Diagonal(ones(N+1))) + alpha * Re * 1im * Diagonal(U_double_prime) * D2

    eigenvalues, eigenvectors = eigen(A)
    return eigenvalues, eigenvectors
end

N = 500   
Re = 2000  # Reynolds number
alpha = 1
c = 0.2
h = 1 / (N + 1)
x = LinRange(0,1-h, N+1)

# Base flow
U = x
U_double_prime = zeros(size(x))
eigenvalues, eigenvectors = orr_sommerfeld_solver(N, Re, alpha, c, U, U_double_prime)

scatter(real(eigenvalues), imag(eigenvalues), label="Eigenvalues", xlabel="Real Part", ylabel="Imaginary Part", legend=:topleft)
title!("Eigenvalues of Orr-Sommerfeld Equation")
