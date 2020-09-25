# Figure 7a: oscillations
using DifferentialEquations, RCall, Roots, SpecialFunctions

# Solution parameters and output
tmax = 400
tspan = (0.0,tmax)
output = Array{Float64}(undef, 0, 4)

# parameters
ι = 1
ℓ = 1
ω = 0.85
φ = 0.09
μ = 0.7
σ² = 0.04
K = 100
q̂ = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/sqrt(2*σ²))) - x, 0, 1))
q = 0.5*(1 + erf((q̂ - μ)/sqrt(2*σ²)))

peq = find_zeros(x -> 0.5*(1 + erf((q + (x-q)*ℓ/(ℓ+ω) - μ)/sqrt(2*σ²))) - x, 0, 1)
yeq = ℓ/(ℓ+ω)

ℛ₀ = ι./(ω*(q̂ .- q .+ yeq*(q .- peq))*yeq)
θ = ℛ₀ .- (1 + sqrt(1 + 2*(1+ι/φ)))/(1 + φ/ι) # threshold for stability
Seq = K ./ ℛ₀ # S equilibrium
Ieq = K*(ℛ₀.-1)./(ℛ₀ .+ ι/φ) # I equilibrium

# ODE
function ODE4comp!(du,u,p,t)
    du[1] = φ*(K-u[1]-u[2]) - ι*u[1]*u[2]/K # Ṡ
    du[2] = ι*u[1]*u[2]/K - ω*(q̂ - q + u[3]*(q-u[4]))*u[2]*u[3] # Ṅ
    du[3] = (1 - u[3])*(q̂ - q + u[3]*(q-u[4]))*(ℓ - ω*u[3]) - ι*u[3]*u[1]/K # ẏ
    du[4] = 0.5*(1 + erf((u[3]*(u[4]-q) + q - μ)/sqrt(2*σ²))) - u[4] # ṗ
end

count = 1
for S₀ = 20:60:80
    I₀ = K - S₀
    for y₀ = 0.2:0.6:0.8
        for p₀ = 0.2:0.6:0.8
            # further initial conditions
            u₀ = [S₀;I₀;y₀;p₀]

            prob = ODEProblem(ODE4comp!,u₀,tspan)
            sol = solve(prob,saveat=0:.01:tmax)
            tsI = sol[2,:]/K
            global output = vcat(output,[collect(0:.01:tmax) tsI count*ones(length(0:.01:tmax),1) zeros(length(0:.01:tmax),1)])
            global count += 1
        end
    end
end

# Initial conditions perturbed about the stable interior equilibrium
ϵ = 0.015
for S₀ = Seq[3]*(1-ϵ):Seq[3]*ϵ:Seq[3]*(1+ϵ) # Seq = 8.5
    for I₀ = Ieq[3]*(1-ϵ):Ieq[3]*ϵ:Ieq[3]*(1+ϵ) # Ieq = 47
        for y₀ = yeq*(1-ϵ):yeq*ϵ:yeq*(1+ϵ) # yeq = 0.54
            for p₀ = peq[3]*(1-ϵ):peq[3]*ϵ:peq[3]*(1+ϵ) # peq = 0.7
                # further initial conditions
                u₀ = [S₀;I₀;y₀;p₀]

                prob = ODEProblem(ODE4comp!,u₀,tspan)
                sol = solve(prob,saveat=0:.01:tmax)
                tsI = sol[2,:]/K
                global output = vcat(output,[collect(0:.01:tmax) tsI count*ones(length(0:.01:tmax),1) ones(length(0:.01:tmax),1)])
                global count += 1
            end
        end
    end
end

@rput output
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

output <- as.data.frame(output)

p <- ggplot(data=output, aes(x=V1, y=V2, group = V3, color = factor(V4))) + geom_line(size=0.2) +
        scale_y_continuous(expand = c(0, 0),breaks=c(0,0.25,0.5,0.75,1),labels=c("0","K/4","K/2","3K/4","K"),lim = c(0,1.01)) +
        scale_x_continuous(expand = c(0, 0),lim = c(0,405)) +
        scale_color_manual(values=c("0"="black","1"="forestgreen")) +
        ylab(expression(paste("Community size, ", I))) +
        xlab("Time") +
        ggtitle("Time series") +
        theme(legend.position="none")

save_plot(p,filename="oscillationsA.png",base_height = 3.5,base_width = 5.25)
"""
