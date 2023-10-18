cd(@__DIR__)

using Distributions
using Printf
using JLD2
using Random
using FFTW
using StaticArrays

# Random.seed!(parse(Int, ARGS[3]))
# CUDA.seed!(parse(Int, ARGS[3]))

# const L = parse(Int, ARGS[1]) # must be a multiple of 4
const L = 16
const λ = 4.0e0
const Γ = 1.0e0
const η = 1.0e0
const T = 1.0e0

const Δt = 0.04e0/Γ
const Rate_phi = Float64(sqrt(2.0*Δt*Γ))
const Rate_pi = Float64(sqrt(2.0*Δt*η))
ξ = Normal(0.0e0, 1.0e0)

function hotstart(n)
	rand(ξ, n, n, n, 3)
end

# nearest neighbor +
function NNp(n)
    n%L+1
end

# nearest neighbor -
function NNm(n)
    (n+L-2)%L+1
end

function pi_step(π, μ, x1, x2)
	norm = cos(2pi*rand())*sqrt(-2*log(rand()))
	q = Rate_pi * norm

    @inbounds δH = 2q * (π[x1..., μ] - π[x2..., μ]) + 2q^2
	P = min(1.0f0, exp(-δH))
	r = rand()

	@inbounds π[x1..., μ] += q * (r<P)
	@inbounds π[x2..., μ] -= q * (r<P)
end

function pi_sweep(π, n, m, μ, (i,j,k))
    xyz = ((2i + m)%L+1, j%L+1, k%L+1)
    @inbounds x1 = (xyz[n%3+1], xyz[(n+1)%3+1], xyz[(n+2)%3+1])
    @inbounds x2 = ((x1[1]-(n!=0))%L+1, (x1[2]-(n!=1))%L+1, (x1[3]-(n!=2))%L+1)

    pi_step(π, μ, x1, x2)
end

function ΔH_phi(x, ϕ, q, m²)
	@inbounds ϕold = ϕ[x...]
	ϕt = ϕold + q
	Δϕ = ϕt - ϕold
	Δϕ² = ϕt^2 - ϕold^2

    @inbounds ∑nn = ϕ[NNp(x[1]), x[2], x[3]] + ϕ[x[1], NNp(x[2]), x[3]] + ϕ[x[1], x[2], NNp(x[3])] + ϕ[NNm(x[1]), x[2], x[3]] + ϕ[x[1], NNm(x[2]), x[3]] + ϕ[x[1], x[2], NNm(x[3])]

	return 3Δϕ² - Δϕ * ∑nn + 0.5m² * Δϕ² + 0.25λ * (ϕt^4 - ϕold^4)
end

function phi_step(m², ϕ, x1, x2)
	norm = cos(2pi*rand())*sqrt(-2*log(rand()))
	q = Rate_phi * norm

    δH = ΔH_phi(x1, ϕ, q, m²) + ΔH_phi(x2, ϕ, -q, m²) + q^2
	P = min(1.0f0, exp(-δH))
	r = rand()

	@inbounds ϕ[x1...] += q * (r<P)
	@inbounds ϕ[x2...] -= q * (r<P)
end

function phi_sweep(m², ϕ, n, m, (i,j,k))
    xyz = ((4i + 2j + m%2)%L+1, (j + k + m÷2)%L+1, k%L+1)
    @inbounds x1 = (xyz[n%3+1], xyz[(n+1)%3+1], xyz[(n+2)%3+1])
    @inbounds x2 = ((x1[1]-(n!=0))%L+1, (x1[2]-(n!=1))%L+1, (x1[3]-(n!=2))%L+1)

    phi_step(m², ϕ, x1, x2)
end

function dissipative(ϕ, π)
    # pi update
    for n in 0:2, m in 0:1
        Threads.@threads for index in 0:3*L^3÷2-1
            μ = index ÷ (L^3÷2) + 1
            i = (index ÷ L^2) % (L ÷ 2)
            j = (index ÷ L) % L
            k = index % L

            pi_sweep(π, n, m, μ, (i,j,k))
        end
    end

    # phi update
    for n in 0:2, m in 0:3
        Threads.@threads for index in 0:L^3÷4-1
            i = index ÷ L^2
            j = (index ÷ L) % L
            k = index % L

            phi_sweep(m², ϕ, n, m, (i,j,k))
        end
    end
end

function Δ(ϕ, n1, n2, n3)
    out = ϕ[NNp(n1), n2, n3] + ϕ[NNm(n1), n2, n3] + ϕ[n1, NNp(n2), n3] + ϕ[n1, NNm(n2), n3] + ϕ[n1, n2, NNp(n3)] + ϕ[n1, n2, NNm(n3)] - 6*ϕ[n1, n2, n3]

    return out
end

function ∇ᵣ(ϕ, n)
    out = ϕ[NNp(n)] - ϕ[n]
    return out 
end

function ∇ₗ(ϕ, n)
    out = - ϕ[NNm(n)] + ϕ[n]
    return out 
end

function Poison(ϕ)
    Afft = fft(ϕ)
    for n1 in 1:L
        for n2 in 1:L
            for n3 in 1:L
                k2 = 4* (sin(pi/L * (n1-1))^2 + sin(pi/L * (n2-1))^2 + sin(pi/L * (n3-1))^2)
                Afft[n1,n2,n3] = -Afft[n1,n2,n3] /  k2 
            end
        end
    end
    Afft[1,1,1] = 0.0 
    ϕ .= real.(ifft(Afft))
end

function projectX(π)
    π_sum = zeros(Float64,(L,L,L))
    temp = zeros(Float64,(L,L,L,3))
    
    for n1 in 1:L, n2 in 1:L, n3 in 1:L
        π_sum[n1,n2,n3] = ∇ᵣ(π[:,n2,n3,1], n1) + ∇ᵣ(π[n1,:,n3,2], n2) + ∇ᵣ(π[n1,n2,:,3], n3)
    end

    for n1 in 1:L, n2 in 1:L, n3 in 1:L 
        temp[n1,n2,n3,1] = - ∇ₗ(@view(π_sum[:,n2,n3]),n1) + Δ(@view(π[:,:,:,1]), n1, n2, n3)

        temp[n1,n2,n3,2] = - ∇ₗ(@view(π_sum[n1,:,n3]),n2) + Δ(@view(π[:,:,:,2]), n1, n2, n3)

        temp[n1,n2,n3,3] = - ∇ₗ(@view(π_sum[n1,n2,:]),n3) + Δ(@view(π[:,:,:,3]), n1, n2, n3)
    end

    for μ in 1:3
        Poison(@view(temp[:,:,:,μ]))    
    end 

    π .= temp

    return 1
end

function thermalize(π, N)
	for _ in 1:N
		dissipative(π)
	end
    projectX(π)
end

m² = -2.28587

π = hotstart(L)

for μ in 1:3
    π[:,:,:,μ] .= π[:,:,:,μ] .- shuffle(π[:,:,:,μ]);
end

@btime projectX(π)

# beforeP= ∇ᵣ(π[:,7,7,1], 7) + ∇ᵣ(π[7,:,7,2], 7) + ∇ᵣ(π[7,7,:,3], 7)
#
# thermalize(π, 1)
#
# afterP1= ∇ᵣ(π[:,7,7,1], 7) + ∇ᵣ(π[7,:,7,2], 7) + ∇ᵣ(π[7,7,:,3], 7)
