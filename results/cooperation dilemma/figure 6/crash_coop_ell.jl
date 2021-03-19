# Figure 6b: equilibrium community size I vs ℓ for the cooperation dilemma
using NLsolve, RCall, Roots, SpecialFunctions

# parameters
ι = 0.1
ω = 1
φ = 1
μ = 0.7 # mean norm sensitivity
σ² = 0.04 # variance in norm sensitivity
K = 100 # total population
q̂ = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/sqrt(2*σ²))) - x, 0, 1))
q = 0.5*(1 + erf((q̂ - μ)/sqrt(2*σ²)))

# output
loweq = Array{Float64}(undef, 0, 3)
mideq = Array{Float64}(undef, 0, 3)
higheq = Array{Float64}(undef, 0, 3)
crasheq = Array{Float64}(undef, 0, 3)

function routhHurwitz(ℓ,p)

    # functions
    f(x) = exp.(-((x .- μ).^2)/(2*σ²))/sqrt(2*π*σ²)

    # equilibria
    y = ℓ/(ℓ+ω)
    p̄ = y*(p .- q) .+ q

    ℛ₀ = ι./(ω*(q̂ .- p̄)*y)
    p = p[ℛ₀ .> 1]
    p̄ = p̄[ℛ₀ .> 1]
    ℛ₀ = ℛ₀[ℛ₀ .> 1]

    S = K./ℛ₀
    I = K*(ℛ₀ .- 1)./(ℛ₀ .+ ι/φ)

    # elements of the Jacobian
    j₁ = φ .+ ι*I/K
    j₂ = φ .+ ι*S/K
    j₃ = ι*I/K
    j₄ = ω*(q̂ .- q .+ 2*(q .- p).*y).*I
    j₅ = ω*I.*y.^2
    j₆ = ι*y/K
    j₇ = (q̂ - q)*(ℓ .+ ω*(1 .- y)) .- (q .- p).*(ℓ .- 2*ℓ*y .- 2*ω.*(1 .- y).*y)
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

# test to see if the population can crash
function crashTest(ℓ)
    # functions
    soly(x) = (μ.-q .+ sqrt(2*σ²)*erfinv.(2*x .- 1))./(x.-q)
    solỹ(x) = -μ.+q .- sqrt(2*σ²)*erfinv.(2*x .- 1) # (q-p)*soly

    # compute crash equilibria y* and p* within [0,1]
    p = find_zeros(p -> ω*solỹ(p)*soly(p)^2 + ω*(q̂-q)*soly(p)^2 - (ℓ+ω)*solỹ(p)*soly(p) - ((ℓ+ω)*(q̂-q) + ι)*soly(p) + ℓ*solỹ(p) + ℓ*(q̂-q), 0, q, no_pts=50)
    p = p[0 .<= p .<= 1]
    y = soly(p)

    p = p[0 .<= y .<= 1]
    y = y[0 .<= y .<= 1]

    ℛ₀ = ι./(ω.*(q̂ .- q .+ y.*(q .- p)).*y)

    # check if the conditions are met
    if any(x->x<1,ℛ₀)
        return true
    else
        return false
    end
end

for m = 1:1:1000
    ℓ = 2e-3*m

    # check if the conditions are met
    if crashTest(ℓ)
        global crasheq = vcat(crasheq, [ℓ 0 0])
    end

    # compute other equilibria
    eqy = ℓ/(ℓ+ω)
    eqp = find_zeros(p -> 0.5*(1 .+ erf.((eqy.*(p.-q) .+ q .- μ)/sqrt(2*σ²))) .- p, 0, 1)
    eqp = eqp[0 .<= eqp .< 1]
    if isempty(eqp)
            break
    end
    eqp = sort(eqp, rev=true)

    # determine p̄ (average p), ℛ₀, and I
    p̄ = eqy.*(eqp .-  q) .+ q
    ℛ₀ = ι./(ω.*(q̂ .- q .+ eqy.*(q .- eqp)).*eqy)
    I = K*(ℛ₀ .- 1)./(ℛ₀ .+ ι./φ)

    # keep values for I>0
    eqp = eqp[I .> 0]
    p̄ = p̄[I .> 0]
    ℛ₀ = ℛ₀[I .> 0]
    I = I[I .> 0]

    if ~isempty(eqp)
        if routhHurwitz(ℓ,eqp[1])
            global higheq = vcat(higheq, [ℓ I[1] 1])
        else
            global higheq = vcat(higheq, [ℓ I[1] 2])
        end
        if length(eqp) > 1
            if routhHurwitz(ℓ,eqp[2])
                global mideq = vcat(mideq, [ℓ I[2] 3])
            else
                global mideq = vcat(mideq, [ℓ I[2] 4])
            end
        end
        if length(eqp) > 3
            if routhHurwitz(ℓ,eqp[3])
                global loweq = vcat(loweq, [ℓ I[3] 5])
            else
                global loweq = vcat(loweq, [ℓ I[3] 6])
            end
        end
    end
end

output = vcat(vcat(loweq,mideq),higheq)
bifpoints = vcat(crasheq[1,1:2]',[crasheq[1,1] higheq[round(Int,0.5*1000*crasheq[1,1]),2]])
bifline = vcat([bifpoints[1,:]' -1],[bifpoints[2,:]' -1])
output = vcat(output,bifline)
output = vcat(output,crasheq)

@rput output bifpoints
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

output <- as.data.frame(output)
bifpoints <- as.data.frame(bifpoints)
colstype = c("-1"="dotted","0"="solid","1"="solid","2"="dashed","3"="solid","4"="dashed","5"="solid","6"="dashed")
colscolour = c("-1"="black","0"="black","1"="black","2"="magenta","3"="black","4"="magenta","5"="solid","6"="dashed")

p <- ggplot() +
        geom_line(data=output, aes(x=V1, y=V2, group = factor(V3), colour = factor(V3), linetype=factor(V3)),size=1) +
        geom_point(data=bifpoints, aes(x=V1, y=V2), shape=21, fill="magenta",color="white") +
        scale_y_continuous(expand = c(0, 0), breaks=c(0,25,50,75,100), labels=c("0","K/4","K/2","3K/4","K"), lim = c(-1.5,101)) +
        scale_x_continuous(expand = c(0, 0),breaks=c(0,0.5,1,1.5,2),labels=c(0,0.5,1,1.5,2),lim = c(0,2.05)) +
        scale_linetype_manual(values=colstype) +
        scale_color_manual(values=colscolour) +
        ylab(expression(paste("Number of insiders, I*"))) +
        xlab(expression(paste("Learning rate, ", '\u2113'))) +
        theme(legend.position = "none") +
        ggtitle(expression(paste("Community size")))

save_plot(p,filename="crash_coop_ell_iota01omega1phi1.png",base_height = 3.5,base_width = 3.5)
"""
