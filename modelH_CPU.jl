#= 

.88b  d88.  .d88b.  d8888b. d88888b db           db   db 
88'YbdP`88 .8P  Y8. 88  `8D 88'     88           88   88 
88  88  88 88    88 88   88 88ooooo 88           88ooo88 
88  88  88 88    88 88   88 88~~~~~ 88           88~~~88 
88  88  88 `8b  d8' 88  .8D 88.     88booo.      88   88 
YP  YP  YP  `Y88P'  Y8888D' Y88888P Y88888P      YP   YP 

=# 

cd(@__DIR__)

using Distributions
using Printf
using Random
using FFTW
using Plots

const LEFT = -1
const RIGHT = 1

const L = 16
const λ = 4.0e0
const Γ = 1.0e0
const η = 1.0e0
const T = 1.0e0
const ρ = 1.0e0

const Δt = 0.04e0/Γ
const Rate_phi = Float64(sqrt(2.0*Δt*Γ))
const Rate_pi = Float64(sqrt(2.0*Δt*η))
ξ = Normal(0.0e0, 1.0e0)

function hotstart(n, Ncomponents)
	rand(ξ, n, n, n, Ncomponents)
end

function hotstart(n)
	rand(ξ, n, n, n)
end

## nearest neighbor
function NNp(n)
    n%L+1
end

function NNm(n)
    (n+L-2)%L+1
end
##

# generates shifts
function shifts(direction, component)
    # quite obscure set of operations to avoid ifs 
    (-direction*(1÷component), -direction*mod(2÷component,2), -direction*(component÷3))
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

function Δ(ϕ)
    # unavoidable allocations
    dϕ = -6*ϕ

    for μ in 1:3
        dϕ .+= circshift(ϕ, shifts(LEFT, μ)) 
        dϕ .+= circshift(ϕ, shifts(RIGHT,μ)) 
    end

    return dϕ
end 

function ∇(src, dest, direction, component) 
    circshift!(dest, src, shifts(direction,component))
    dest .= (dest .- src) * direction

    return 1 
end

function ∇conj(src, dest, direction, component) 
    circshift!(dest, src, shifts(-direction,component))
    dest .= (src .- dest) * direction

    return 1 
end

function Poisson(ϕ)
    Afft = fft(ϕ)
    for n1 in 1:L, n2 in 1:L, n3 in 1:L
        k2 = 4* (sin(pi/L * (n1-1))^2 + sin(pi/L * (n2-1))^2 + sin(pi/L * (n3-1))^2)
        Afft[n1,n2,n3] = -Afft[n1,n2,n3] /  k2 
    end
    Afft[1,1,1] = 0.0 
    ϕ .= real.(ifft(Afft))
end

function project(π, direction)
    # ∇_μ projectμ = 0 

    π_sum = zeros(Float64,(L,L,L))
    temp = zeros(Float64,(L,L,L,3))

    for μ in 1:3 
        ∇(@view(π[:,:,:,μ]), @view(temp[:,:,:,μ]), direction, μ)
        π_sum .-= temp[:,:,:,μ] # note that π_sum is actiually - π_sum
    end 

    Threads.@threads for μ in 1:3 
        ∇conj(π_sum, @view(temp[:,:,:,μ]), direction, μ)
        temp[:,:,:,μ] .+= Δ(@view(π[:,:,:,μ]))
        Poisson(@view(temp[:,:,:,μ]))    
    end 

    π .= temp

    return 1
end

# maccormack
function maccormack(ϕ, π, direction)
    # project(π, direction)

    # temporary arrays 
    dj = similar(ϕ)
    dϕ = similar(ϕ)
    dϕ_μ = similar(ϕ)
    dϕ_ν = similar(ϕ)
    dπ = similar(π)

    dj .= 0.0
    dϕ .= ϕ
    dπ .= 0.0

    for μ in 1:3
        ∇(π[:,:,:,μ] .* ϕ, dj, direction, μ)
        dϕ .-= Δt*dj        
    end

    dj .= 0.0

    for ν in 1:3
        ∇conj(ϕ, dϕ_ν, direction, ν)

        for μ in 1:3
            ∇conj(ϕ, dϕ_μ, direction, μ) 
            dj .= 1/ρ * π[:,:,:,μ] .* π[:,:,:,ν] .+ dϕ_μ .* dϕ_ν
            # overwrite dϕ_μ because it's not needed
            ∇(dj, dϕ_μ, direction, μ)
            dπ[:,:,:,ν] .-= Δt*dϕ_μ 
        end

        dπ[:,:,:,ν] .+= π[:,:,:,ν] 
    end

    return (dϕ,dπ) # returns * or ** updates 
end

function maccormack_step(ϕ, π)
    step1 = maccormack(ϕ, π, RIGHT)
    step2 = maccormack(step1..., LEFT)
    ϕ .= 1/2 * (ϕ .+ step2[1])
    π .= 1/2 * (π .+ step2[2])
end

function thermalize(ϕ, π, N)
	for _ in 1:N
		dissipative(ϕ, π)
	end
end

m² = -2.28587

π = hotstart(L, 3)
ϕ = hotstart(L)

# for μ in 1:3
#     π[:,:,:,μ] .= π[:,:,:,μ] .- shuffle(π[:,:,:,μ]);
# end

π .= 0.0e0
π[:,:,:,1] .= 0.1e0
ϕ .= [exp(-((x-8)^2 + (y-8)^2 + (z-8)^2) / 2) for x in 1:L, y in 1:L, z in 1:L]

# plot(π[:,8,8,1])
plot(ϕ[:,8,8])
println(sum(ϕ))

for _ in 1:70
    maccormack_step(ϕ, π)
end

println(sum(ϕ))
plot!(ϕ[:,8,8])
