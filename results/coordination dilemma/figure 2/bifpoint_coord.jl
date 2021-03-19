# Figure 2b: bifurcation points for various σ²
using NLsolve, RCall, Roots, SpecialFunctions

μ = 0.5 # mean norm sensitivity
output = zeros(1000,3)

for m = 1:1:1000
    σ² = 1e-4*m # variance in norm sensitivity

    q̂ = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/sqrt(2*σ²))) - x, 0, 1))
    q = 0.5*(1+erf((q̂-μ)/sqrt(2*σ²)))

    # compute roots where x[1] = p and x[2] = τ = ω/ℓ
    function g!(G,x)
        G[1] = 0.5*(1 + erf(((x[1]-q)/(1+x[2]) + q - μ)/sqrt(2*σ²))) - x[1]
        G[2] = exp(-(((x[1]-q)/(1+x[2]) + q - μ)^2)/(2*σ²))/sqrt(2*π*σ²)-1-x[2]
    end

    function j!(J,x)
        J[1,1] = exp(-(((x[1]-q)/(1+x[2]) + q - μ)^2)/(2*σ²))/(sqrt(2*π*σ²)*(1+x[2])) - 1
        J[1,2] = -(x[1]-q)*exp(-(((x[1]-q)/(1+x[2]) + q - μ)^2)/(2*σ²))/(sqrt(2*π*σ²)*(1+x[2])^2)
        u = exp(-(((x[1]-q)/(1+x[2]) + q - μ)^2)/(2*σ²))/sqrt(2*π*σ²)
        J[2,1] = (1/(1+x[2])^2)*(-((x[1]-q)/(1+x[2]) + q - μ)/σ²)*u
        J[2,2] = -((x[1]-q)/(1+x[2])^3)*(-((x[1]-q)/(1+x[2]) + q - μ)/σ²)*u
    end

    # initial guesses for root finder
    p₀ = 0.4
    τ₀ = 0.1 # bifurcation point, τ = ω/ℓ

    sol = nlsolve(g!, j!, [p₀; τ₀]) # find roots

    global output[m,:] = [σ² sol.zero[2] sol.zero[1]]
end

output = output[output[:,2] .< 1,:]

@rput output
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

output <- as.data.frame(output)

p <- ggplot(data=output, aes(x=V1, y=V2)) + geom_line(size=1) +
        scale_y_continuous(expand = c(0, 0),breaks=c(0,0.25,0.5,0.75,1),labels=c(0,0.25,0.5,0.75,1),lim = c(0,1.01)) +
        scale_x_continuous(expand = c(0, 0),labels=c(0,0.025,0.05,0.075,0.1),lim = c(0,0.105)) +
        ylab(expression(paste("Ratio of outflow to learning rates, ", omega,'/\u2113'))) +
        xlab(expression(paste("Variance in norm sensitivity, ",sigma^2))) +
        ggtitle(expression("Bifurcation points")) +
        theme(axis.title.x=element_text(margin=margin(t=0,r=0,b=-0.75,l=0)))

save_plot(p,filename="bifpoint_coord.png",base_height = 3.5,base_width = 3.5)
"""
