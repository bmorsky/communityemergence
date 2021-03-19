# Figure 4c: equilibria ṗ=0 vs bifurcation parameter ω/ℓ
using RCall, Roots, SpecialFunctions

# parameters
μ = 0.7 # mean norm sensitivity
σ² = 0.09 # variance in norm sensitivity
K = 100 # total population
q̂ = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/sqrt(2*σ²))) - x, 0, 1))
q = 0.5*(1 + erf((q̂ - μ)/sqrt(2*σ²)))

loweq = Array{Float64}(undef, 0, 3)
mideq = Array{Float64}(undef, 0, 3)
higheq = Array{Float64}(undef, 0, 3)

for m = 1:1:10000
    τ = 2e-4*m # bifurcation parameter, τ = ω/ℓ

    # compute p equilibrium
    eqy = 1/(1+τ)
    eqp = find_zeros(p -> 0.5*(1 .+ erf.((eqy.*(p.-q) .+ q .- μ)/sqrt(2*σ²))) .- p, 0, 1)
    eqp = eqp[0 .<= eqp .< 1]
    if isempty(eqp)
            break
    end
    p̄ = eqy.*(eqp .-  q) .+ q
    global higheq = vcat(higheq, [τ maximum(p̄) 3])
    if length(eqp) > 1
            global loweq = vcat(loweq, [τ p̄[1] 1])
            global mideq = vcat(mideq, [τ p̄[2] 2])
    end
end

output = higheq

@rput output
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

output <- as.data.frame(output)

p <- ggplot() +
        geom_line(data=output, aes(x=V1, y=V2, group = factor(V3),linetype=factor(V3)),size=1) +
        scale_y_continuous(expand = c(0, 0),breaks=c(0,0.5,1),labels=c(0,0.5,1),lim = c(0,1.01)) +
        scale_x_continuous(expand = c(0, 0),breaks=c(0,0.5,1,1.5,2),labels=c(0,0.5,1,1.5,2),lim = c(0,2.01)) +
        scale_linetype_manual(values=c("solid","solid","solid")) +
        ylab(expression(paste("Mean frequency of cooperation, ", bar(p),"*"))) +
        xlab(expression(paste("Ratio of outflow to learning rates, ", omega,'/\u2113'))) +
        theme(legend.position = "none", axis.title.x = element_text(hjust=0.8)) +
        ggtitle(expression(paste("Equilibria ",dot(p),"=0 with ",sigma^2,"=0.09")))

save_plot(p,filename="pEq_coop_var09.png",base_height = 3.5,base_width = 3.5)
"""
