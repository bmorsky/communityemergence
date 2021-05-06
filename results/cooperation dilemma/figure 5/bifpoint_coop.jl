# Figure 4d: bifurcation points for various σ²
using RCall, Roots, SpecialFunctions

# parameters
μ = 0.7 # mean norm sensitivity
K = 100 # total population

# output
output1 = zeros(1000,3)
output2 = zeros(100,3)

for m = 1:1:1000
        σ² = 1e-5 + 4.2e-5*m # variance in norm sensitivity
        q̂ = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/sqrt(2*σ²))) - x, 0, 1))
        q = 0.5*(1 + erf((q̂ - μ)/sqrt(2*σ²)))
        eq = Array{Float64}(undef, 0, 2)
        for p̄ = 0.3:0.01:0.75
                # calculate the equilibrium threshold τ = ω/ℓ
                eqτ = find_zeros(τ -> 0.5*(1 .+ erf.(((p̄ .+ (p̄.-q).*τ.-q)./(1+τ) .+ q .- μ)/sqrt(2*σ²))) .- p̄ .- (p̄.-q).*τ, 0, 1000)
                eqτ = eqτ[0 .<= eqτ]
                if isempty(eqτ)
                        break
                end
                eq = vcat(eq, [p̄ eqτ])
        end
        if ~isempty(eq)
                global output1[m,:] = [σ² maximum(eq[:,2]) 0]
        end
end

for m = 1:1:100
        σ² = 0.018*(m-1)/100 + 0.0244
        q̂ = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/(sqrt(2*σ²)))) - x, 0, 1))
        q = 0.5*(1 + erf((q̂ - μ)/(sqrt(2*σ²))))
        eq = Array{Float64}(undef, 0, 2)
        for p̄ = 0.6:0.01:0.9
                # calculate the equilibrium threshold τ = ω/ℓ
                eqτ = find_zeros(τ -> 0.5*(1 .+ erf.(((p̄ .+ (p̄.-q).*τ.-q)./(1+τ) .+ q .- μ)/(sqrt(2*σ²)))) .- p̄ .- (p̄.-q).*τ, 0, 100)
                eqτ = eqτ[0 .<= eqτ]
                if isempty(eqτ)
                        break
                end
                eq = vcat(eq, [p̄ eqτ])
        end
        global output2[m,:] = [σ² minimum(eq[:,2]) 1]
end

output = vcat(output1,output2)

@rput output
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

output <- as.data.frame(output)

p <- ggplot(data=output, aes(x=V1, y=V2, group = factor(V3))) + geom_line(size=1) +
        scale_y_continuous(expand = c(0, 0),breaks=c(0,0.5,1,1.5,2),labels=c(0,0.5,1,1.5,2),lim = c(0,2.02)) +
        scale_x_continuous(expand = c(0, 0),breaks=c(0,0.025,0.05,0.075,0.1),labels=c(0,0.025,0.05,0.075,0.1),lim = c(0,0.1)) +
        ylab(expression(paste("Ratio of outflow to learning rates, ", omega,'/\u2113'))) +
        xlab(expression(paste("Variance in norm sensitivity, ",sigma^2))) +
        theme(legend.position = "none") +
        ggtitle(expression("Bifurcation points")) +
        theme(axis.title.x=element_text(margin=margin(t=0,r=0,b=-0.5,l=0)))

save_plot(p,filename="bifpoint_coop.png",base_height = 3.5,base_width = 3.5)
"""
