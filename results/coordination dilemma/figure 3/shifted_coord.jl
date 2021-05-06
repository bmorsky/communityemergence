# Figure 2c: shifted norm sensitivity for the coordination dilemma
using NLsolve, RCall, Roots, SpecialFunctions

p = collect(0:0.001:1)
μ₅ = 0.5 # mean norm sensitivity
σ² = 0.04 # variance in norm sensitivity

function bifpoint(q,μ,σ²)
    # compute roots
    function g!(G,x)
        G[1] = 0.5*(1 + erf(((x[1]-q)/(1+x[2]) + q - μ)/sqrt(2*σ²))) - x[1]
        G[2] = exp(-(((x[1]-q)/(1+x[2]) + q - μ)^2)/(2*σ²))/(sqrt(2*π*σ²)*(1+x[2])) - 1
    end

    function j!(J,x)
        J[1,1] = exp(-(((x[1]-q)/(1+x[2]) + q - μ)^2)/(2*σ²))/(sqrt(2*π*σ²)*(1+x[2])) - 1
        J[1,2] = -(x[1]-q)*exp(-(((x[1]-q)/(1+x[2]) + q - μ)^2)/(2*σ²))/(sqrt(2*π*σ²)*(1+x[2])^2)
        u = exp(-(((x[1]-q)/(1+x[2]) + q - μ)^2)/(2*σ²))/(sqrt(2*π*σ²)*(1+x[2]))
        J[2,1] = (1/(1+x[2]))*(-((x[1]-q)/(1+x[2]) + q - μ)/σ²)*u
        J[2,2] = -((x[1]-q)/(1+x[2])^2)*(-((x[1]-q)/(1+x[2]) + q - μ)/σ²)*u
    end

    return nlsolve(g!, j!, [0.1; 0.1]).zero
end

q₅ = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/sqrt(2*σ²))) - x, 0, 1))

(p₅, τ₅) = bifpoint(q₅,μ₅,σ²) # find bifurcation point

F₅ = 0.5*(1 .+ erf.((p .- μ₅)/sqrt(2*σ²)))
F₅shift = 0.5*(1 .+ erf.(((p .- q₅)./(1 .+ τ₅) .+ q₅ .- 0.5)/sqrt(2*σ²)))

output = vcat([p p zeros(length(p),1)],vcat([p F₅ ones(length(p),1)],[p F₅shift 2*ones(length(p),1)]))

@rput output
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

output <- as.data.frame(output)

p <- ggplot(data=output, aes(x=V1, y=V2, group = factor(V3), color=factor(V3))) + geom_line(aes(linetype=factor(V3)),size=1) +
        scale_y_continuous(expand = c(0, 0),breaks=c(0,0.25,0.5,0.75,1),labels=c(0,0.25,0.5,0.75,1),lim = c(0,1.01)) +
        scale_x_continuous(expand = c(0, 0),labels=c(0,0.25,0.5,0.75,1),lim = c(0,1.05)) +
        scale_linetype_manual(values=c("solid","solid","dashed")) +
        scale_color_manual(values=c("black","blue","blue")) +
        ylab(expression(paste("Conditional cooperation, ",  F(bar(p))))) +
        theme(legend.position = "none") +
        xlab(expression(paste("Frequency of cooperation, ", bar(p)))) +
        ggtitle(expression("Norm sensitivity"))

save_plot(p,filename="shifted_coord.png",base_height = 3.5,base_width = 3.5)
"""
