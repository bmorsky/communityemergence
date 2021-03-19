# Figure 2a: equilibria ṗ=0 vs bifurcation parameter ω/ℓ
using RCall, Roots, SpecialFunctions

# parameters
μ = 0.5 # mean norm sensitivity
σ² = 0.04 # variance in norm sensitivity
K = 100 # total population
q = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/sqrt(2*σ²))) - x, 0, 1))
f(y,p) = exp.(-(((p.-q).*y .+ q .- μ).^2)/(2*σ²))/sqrt(2*π*σ²)

loweq = Array{Float64}(undef, 0, 3)
mideq = Array{Float64}(undef, 0, 3)

for m = 1:1:1000
    τ = 4e-4*m # bifurcation parameter, τ = ω/ℓ

    # compute p equilibrium
    eqy = 1/(1+τ) # y equilibrium
    eqp = find_zeros(p -> 0.5*(1 .+ erf.((eqy.*(p.-q) .+ q .- μ)/sqrt(2*σ²))) .- p, 0, q) # p equilibrium
    eqp = eqp[0 .<= eqp .< q]
    if isempty(eqp)
            break
    end
    p̄ = eqy.*(eqp .-  q) .+ q
    if eqy.*f(eqy,eqp)[1] < 1
            global loweq = vcat(loweq, [τ minimum(p̄) 1])
    else
            global loweq = vcat(loweq, [τ minimum(p̄) 2])
    end
    if eqy.*f(eqy,eqp)[2] < 1
            global mideq = vcat(mideq, [τ maximum(p̄) 3])
    else
            global mideq = vcat(mideq, [τ maximum(p̄) 4])
    end
end

output = vcat(vcat(loweq,mideq),[collect(1:1:1000)*4e-4 q*ones(1000,1) zeros(1000,1)])
bifpoint = 0.5*(loweq[end,1:2]+mideq[end,1:2])'
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
        geom_line(data=output, aes(x=V1, y=V2, color=factor(V3),group = factor(V3),linetype=factor(V3)),size=1) +
        geom_point(data=bifpoint, aes(x=V1, y=V2), shape=21, fill="magenta",color="white") +
        scale_y_continuous(expand = c(0, 0),breaks=c(0,0.25,0.5,0.75,1),labels=c(0,0.25,0.5,0.75,1),lim = c(0,1.01)) +
        scale_x_continuous(expand = c(0, 0),breaks=c(0,0.1,0.2,0.3,0.4),labels=c(0,0.1,0.2,0.3,0.4),lim = c(0,0.42)) +
        scale_linetype_manual(values=colstype) +
        scale_color_manual(values=colscolour) +
        ylab(expression(paste("Mean frequency of cooperation, ", bar(p),"*"))) +
        xlab(expression(paste("Ratio of outflow to learning rates, ", omega,'/\u2113'))) +
        theme(legend.position = "none", axis.title.x = element_text(hjust=0.8)) +
        ggtitle(expression(paste("Equilibria ",dot(p),"=0")))

save_plot(p,filename="pEq_coord.png",base_height = 3.5,base_width = 3.5)
"""
