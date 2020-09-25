# Figure 3a: equilibria I vs ℓ
using RCall, Roots, SpecialFunctions

# parameters
ι = 1
ω = 0.2
φ = 1
μ = 0.5 # mean norm sensitivity
σ² = 0.04 # variance in norm sensitivity
K = 100 # total population
q = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/sqrt(2*σ²))) - x, 0, 1))

loweq = Array{Float64}(undef, 0, 3)
mideq = Array{Float64}(undef, 0, 3)
higheq = Array{Float64}(undef, 0, 3)

function routhHurwitz(ℓ,p)

    # functions
    f(x) = exp.(-((x .- μ).^2)/(2*σ²))/sqrt(2*π*σ²)

    # equilibria
    y = ℓ/(ℓ+ω)
    p̄ = y*(p .- q) .+ q

    ℛ₀ = ι./(ω*(q .- p)*y.^2)
    p = p[ℛ₀ .> 1]
    p̄ = p̄[ℛ₀ .> 1]
    ℛ₀ = ℛ₀[ℛ₀ .> 1]

    S = K./ℛ₀
    I = K*(ℛ₀ .- 1)./(ℛ₀ .+ ι/φ)

    # elements of the Jacobian
    j₁ = φ .+ ι*I/K
    j₂ = φ .+ ι*S/K
    j₃ = ι*I/K
    j₄ = 2*ω*(q .- p).*y.*I
    j₅ = ω*I.*y.^2
    j₆ = ι*y/K
    j₇ = (p .- q).*(ℓ .- 2*ℓ*y .- 2*ω.*(1 .- y).*y)
    j₈ = y.*(1 .- y).*(ℓ .- ω*y)
    j₉ = (q .- p).*f(p̄)
    j₁₀ = 1 .- y.*f(p̄)

    a₃ = j₁ .+ j₇ .+ j₁₀
    a₂ = j₂.*j₃ .+ j₁.*j₇ .+ j₁.*j₁₀ .+ j₇.*j₁₀ .- j₈.*j₉
    a₁ = j₂.*(j₄.*j₆ .+ j₃.*(j₇ .+ j₁₀)) .+ j₁.*(j₇.*j₁₀ .- j₈.*j₉)
    a₀ = j₂.*(j₆.*(j₅.*j₉ .+ j₄.*j₁₀) .+ j₃.*(j₇.*j₁₀ .- j₈.*j₉))

    Δ₄ = a₁.*(a₂.*a₃ .- a₀) .- a₀.*a₃.^2 .> 0
    Δ₃ = a₃ .> 0
    Δ₂ = a₂ .> 0
    Δ₁ = a₁ .> 0
    Δ₀ = a₀ .> 0

    if Δ₄.*Δ₃.*Δ₂.*Δ₁.*Δ₀ .== 1
        return true
    else
        return false
    end
end

for m = 1:1:1000
    ℓ = 2e-3*m

    # compute p equilibrium
    eqy = ℓ/(ℓ+ω) # y equilibrium
    eqp = find_zeros(p -> 0.5*(1 .+ erf.((eqy.*(p.-q) .+ q .- μ)/sqrt(2*σ²))) .- p, 0, 1) # p equilibrium
    eqp = eqp[0 .<= eqp .< 1]
    if isempty(eqp)
            break
    end
    eqp = sort(eqp, rev=true)
    p̄ = eqy.*(eqp .-  q) .+ q
    ℛ₀ = ι./(ω.*(q .- eqp).*eqy^2)
    I = K*(ℛ₀ .- 1)./(ℛ₀ .+ ι./φ)
    I[findmax(p̄)[2]] = K
    global higheq = vcat(higheq, [ℓ I[1] 0])
    if length(eqp) > 1
            if routhHurwitz(ℓ,eqp[2])
                    global mideq = vcat(mideq, [ℓ I[2] 1])
            else
                    global mideq = vcat(mideq, [ℓ I[2] 2])
            end
            if routhHurwitz(ℓ,eqp[3])
                    global loweq = vcat(loweq, [ℓ I[3] 3])
            else
                    global loweq = vcat(loweq, [ℓ I[3] 4])
            end
    end
end

output = vcat(vcat(loweq,mideq),higheq)
bifpoint = 0.5*(loweq[1,1:2]+mideq[1,1:2])'
bifline = vcat([bifpoint -1],[bifpoint[1] 0 -1])
output = vcat(output,bifline)

@rput output bifpoint
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

output <- as.data.frame(output)
bifpoint <- as.data.frame(bifpoint)
colstype = c("-1"="dotted","0"="solid","1"="solid","2"="dashed","3"="solid","4"="dashed")
colscolour = c("-1"="black","0"="black","1"="black","2"="magenta","3"="black","4"="magenta")

p <- ggplot() +
        geom_line(data=output, aes(x=V1, y=V2, color = factor(V3) ,group = factor(V3),linetype=factor(V3)),size=1) +
        geom_point(data=bifpoint, aes(x=V1, y=V2), shape=21, fill="magenta",color="white") +
        geom_text(data=bifpoint,aes(x=V1, y=V2),label = expression(paste('\u2113',"*")),hjust=1,vjust=0.07) +
        scale_y_continuous(expand = c(0, 0),breaks=c(0,25,50,75,100),labels=c("0","K/4","K/2","3K/4","K"),lim = c(0,101)) +
        scale_x_continuous(expand = c(0, 0),breaks=c(0,0.5,1,1.5,2),labels=c(0,0.5,1,1.5,2),lim = c(0,2.05)) +
        scale_linetype_manual(values=colstype) +
        scale_color_manual(values=colscolour) +
        ylab(expression(paste("Size of community at equilibrium, I*"))) +
        xlab(expression(paste("Learning rate, ", '\u2113'))) +
        theme(legend.position = "none") +
        ggtitle(expression(paste("Community size")))

save_plot(p,filename="IEq_coord_ell_iota1omega02phi1.png",base_height = 3.5,base_width = 3.5)
"""
