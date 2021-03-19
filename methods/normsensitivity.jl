# Figure 1: norm sensitivity for different means
using  NLsolve, RCall, Roots, SpecialFunctions

p = collect(0:0.001:1)
μ₅ = 0.5 # mean norm sensitivity for the coordination dilemma
μ₇ = 0.7 # mean norm sensitivity for the cooperation dilemma
σ² = 0.04 # variance in norm sensitivity
F₅ = 0.5*(1 .+ erf.((p .- μ₅)/sqrt(2*σ²))) # CDF for the coordination dilemma
F₇ = 0.5*(1 .+ erf.((p .- μ₇)/sqrt(2*σ²))) # CDF for the cooperation dilemma

output = vcat([p p zeros(length(p),1)],vcat([p F₅ ones(length(p),1)],[p F₇ 2*ones(length(p),1)]))

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
        scale_linetype_manual(values=c("solid","solid","solid")) +
        scale_color_manual(values=c("black","blue","red")) +
        ylab(expression(paste("Conditional cooperation, ", F(p)))) +
        theme(legend.position = "none") +
        xlab(expression(paste("Frequency of cooperation, ", p))) +
        ggtitle(expression("Norm sensitivity"))

save_plot(p,filename="normsensitivity.png",base_height = 3.5,base_width = 3.5)
"""
