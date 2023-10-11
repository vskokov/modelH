cd(@__DIR__)

using Distributions
using Printf
using JLD2
using Random
using FFTW

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

function project(π)
    π_fft = cat([fft(π[:,:,:,i]) for i in 1:3]...; dims=4)
    
    π_sum = [2*sin(pi/L * (i-1))*π[i,j,k,1] + 2*sin(pi/L * (j-1))*π[i,j,k,2] + 2*sin(pi/L * (k-1))*π[i,j,k,3] for i in 1:L, j in 1:L, k in 1:L]
    
    π_fft[:,:,:,1] .-= [2*sin(pi/L * i) / 4 / (sin(pi/L * i)^2 + sin(pi/L * j)^2 + sin(pi/L * k)^2) for i in 0:L-1, j in 0:L-1, k in 0:L-1] .* π_sum
    π_fft[:,:,:,2] .-= [2*sin(pi/L * j) / 4 / (sin(pi/L * i)^2 + sin(pi/L * j)^2 + sin(pi/L * k)^2) for i in 0:L-1, j in 0:L-1, k in 0:L-1] .* π_sum
    π_fft[:,:,:,3] .-= [2*sin(pi/L * k) / 4 / (sin(pi/L * i)^2 + sin(pi/L * j)^2 + sin(pi/L * k)^2) for i in 0:L-1, j in 0:L-1, k in 0:L-1] .* π_sum
    π_fft[1,1,1,:] .= 0

    cat([ifft(π_fft[:,:,:,i]) for i in 1:3]...; dims=4)
end

function dissipative(π)
    for n in 0:2, m in 0:1
        Threads.@threads for index in 0:3*L^3÷2-1
            μ = index ÷ (L^3÷2)
            i = (index ÷ L^2) % (L ÷ 2)
            j = (index ÷ L) % L
            k = index % L
            
            pi_sweep(π, n, m, μ, (i,j,k))
        end
    end
end

function thermalize(π, N)
	for i in 1:N
		dissipative(π)
	end
end

m² = -2.28587

π = hotstart(L)
π[:,:,:,1] .= π[:,:,:,1] .- shuffle(π[:,:,:,1]);
π[:,:,:,2] .= π[:,:,:,2] .- shuffle(π[:,:,:,2]);
π[:,:,:,3] .= π[:,:,:,3] .- shuffle(π[:,:,:,3]);
