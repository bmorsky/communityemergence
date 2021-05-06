# Figure 4b: equilibria ṗ=0 vs bifurcation parameter ω/ℓ
using RCall, Roots, SpecialFunctions

# parameters
μ = 0.7 # mean norm sensitivity
σ² = 0.04 # variance in norm sensitivity
K = 100 # total population
q̂ = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/sqrt(2*σ²))) - x, 0, 1))
q = 0.5*(1 + erf((q̂ - μ)/sqrt(2*σ²)))

# functions
f(x) = exp.(-((x .- μ).^2)/(2*σ²))/sqrt(2*π*σ²)

eq = Array{Float64}(undef, 0, 3)

counter = 1
for m = 1:1:2000
    p̄ = 5e-3*m

    eqτ = find_zeros(τ -> 0.5*(1 .+ erf.(((p̄ .+ (p̄.-q).*τ.-q)./(1+τ) .+ q .- μ)/sqrt(2*σ²))) .- p̄ .- (p̄.-q).*τ, 0, 1000)
    eqτ = eqτ[0 .<= eqτ]
    if isempty(eqτ)
            break
    end
    if (1/(1+eqτ[1]))*f(p̄) < 1
            global eq = vcat(eq, [p̄ eqτ counter])
    else
            global eq = vcat(eq, [p̄ eqτ 2])
            global counter = 3
    end
end

output = eq
bif1 = findmax(eq[1:150,2])[2]
bif2 = findmin(eq[150:end,2])[2]+150
bifpoints = [eq[bif1,1:2]'; eq[bif2,1:2]']
bifline = vcat([bifpoints [-2;-1]],[[0;0] bifpoints[:,2] [-2;-1]])
output = vcat(output,bifline)

@rput output bifpoints
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

output <- as.data.frame(output)
bifpoints <- as.data.frame(bifpoints)
colstype = c("-2"="dotted","-1"="dotted","0"="solid","1"="solid","2"="dashed","3"="solid")
colscolour = c("-2"="black","-1"="black","0"="black","1"="black","2"="magenta","3"="black")

p <- ggplot() +
        geom_line(data=output, aes(x=V1, y=V2, color = factor(V3), group = factor(V3),linetype=factor(V3)),size=1) +
        geom_point(data=bifpoints, aes(x=V1, y=V2), shape=21, fill="magenta",color="white") +
        scale_y_continuous(expand = c(0, 0),labels=c(0,0.5,1,1.5),lim = c(0,1.55)) +
        scale_x_continuous(expand = c(0, 0),breaks=c(0,0.5,1),labels=c(0,0.5,1),lim = c(0,1.01)) +
        scale_linetype_manual(values=colstype) +
        scale_color_manual(values=colscolour) +
        xlab(expression(paste("Mean frequency of cooperation, ", bar(p),"*"))) +
        ylab(expression(paste("Ratio of outflow to learning rates, ", omega,'/\u2113'))) +
        theme(legend.position = "none", axis.title.x = element_text(hjust=0.8)) +
        ggtitle(expression(paste("Equilibria ",dot(p),"=0 with ",sigma^2,"=0.04"))) +
        coord_flip()

save_plot(p,filename="pEq_coop_var04.png",base_height = 3.5,base_width = 3.5)
"""
