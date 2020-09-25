using RCall, Roots, SpecialFunctions

# parameters
K = 100 # total population
μ = 0.7
σ² = 0.04
q̂ = maximum(find_zeros(x -> 0.5*(1 + erf((x - 0.5)/sqrt(2*σ²))) - x, 0, 1))
q = 0.5*(1 + erf((q̂ - μ)/sqrt(2*σ²)))

function routhHurwitz(μ,σ²,ι,ℓ,ω,φ)

    # functions
    f(x) = exp.(-((x .- μ).^2)/(2*σ²))/sqrt(2*π*σ²)

    # equilibria
    y = ℓ/(ℓ+ω)
    p = find_zeros(x -> 0.5*(1 + erf((y*(x - q) + q - μ)/sqrt(2*σ²))) - x, 0, q)
    p = p[0 .<= p .< 1]
    p = sort(p, rev=true)
    p̄ = y*(p .- q) .+ q

    ℛ₀ = ι./(ω*(q̂ .- p̄)*y)
    p = p[ℛ₀ .> 1]
    p̄ = p̄[ℛ₀ .> 1]
    ℛ₀ = ℛ₀[ℛ₀ .> 1]

    S = K./ℛ₀
    I = K*(ℛ₀ .- 1)./(ℛ₀ .+ ι/φ)
    #i=1
    #eigvals([-j₁[i] -j₂[i] 0 0; j₃[i] 0 -j₄[i] j₅[i]; -j₆ 0 -j₇[i] -j₈; 0 0 -j₉[i] -j₁₀[i]])

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

    return [length(p), length(p[Δ₄.*Δ₃.*Δ₂.*Δ₁.*Δ₀ .== 1])]
end

function crashTest(μ,σ²,ι,ℓ,ω,φ)
    # functions
    soly(x) = (μ.-q .+ sqrt(2*σ²)*erfinv.(2*x .- 1))./(x.-q)
    solỹ(x) = -μ.+q .- sqrt(2*σ²)*erfinv.(2*x .- 1) # (q-p)*soly

    # compute crash equilibria y* and p* within [0,1]
    p = find_zeros(p -> ω*solỹ(p)*soly(p)^2 + ω*(q̂-q)*soly(p)^2 - (ℓ+ω)*solỹ(p)*soly(p) - ((ℓ+ω)*(q̂-q) + ι)*soly(p) + ℓ*solỹ(p) + ℓ*(q̂-q), 0, q, no_pts=350)
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

function categorize(a,b)
    if ~b
        return  a[2]
    else
        return  10+a[2]
    end
end

num = 200
ιVℓ = Array{Float64}(undef, num^2, 3)
ιVω = Array{Float64}(undef, num^2, 3)
ιVφ = Array{Float64}(undef, num^2, 3)

ℓVω = Array{Float64}(undef, num^2, 3)
φVℓ = Array{Float64}(undef, num^2, 3)

φVω = Array{Float64}(undef, num^2, 3)

counter=1
for m = 1:1:num
    for n = 1:1:num
        M₁ = 0.6*m/num
        M₂ = 0.6 + 0.8*m/num
        N₁ = 0.6*n/num
        N₂ = 0.6 + 0.8*n/num
        ιVℓ[counter,:] = [M₁ N₂ categorize(routhHurwitz(μ,σ²,M₁,N₂,1,0.5), crashTest(μ,σ²,M₁,N₂,1,0.5))]
        ιVω[counter,:] = [M₁ N₂ categorize(routhHurwitz(μ,σ²,M₁,1,N₂,0.5), crashTest(μ,σ²,M₁,1,N₂,0.5))]
        ιVφ[counter,:] = [M₁ N₁ categorize(routhHurwitz(μ,σ²,M₁,1,1,N₁), crashTest(μ,σ²,M₁,1,1,N₁))]
        ℓVω[counter,:] = [M₂ N₂ categorize(routhHurwitz(μ,σ²,1,M₂,N₂,0.5), crashTest(μ,σ²,1,M₂,N₂,0.5))]
        φVℓ[counter,:] = [M₁ N₂ categorize(routhHurwitz(μ,σ²,1,N₂,1,M₁), crashTest(μ,σ²,1,N₂,1,M₁))]
        φVω[counter,:] = [M₁ N₂ categorize(routhHurwitz(μ,σ²,1,1,N₂,M₁), crashTest(μ,σ²,1,1,N₂,M₁))]
        global counter = counter + 1
    end
end

iotaVell = ιVℓ
iotaVomega = ιVω
iotaVphi = ιVφ
ellVomega = ℓVω
phiVell = φVℓ
phiVomega = φVω

@rput iotaVell iotaVomega iotaVphi ellVomega phiVell phiVomega
R"""
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(viridis)
library(viridisLite)

iotaVell <- as.data.frame(iotaVell)
iotaVomega <- as.data.frame(iotaVomega)
iotaVphi <- as.data.frame(iotaVphi)
ellVomega <- as.data.frame(ellVomega)
phiVell <- as.data.frame(phiVell)
phiVomega <- as.data.frame(phiVomega)

cols = c("0"="#2EC4b6", "1"="#F5F5F5", "2"="#011627", "10"="#E71D36", "11"="#FF9F1C")
iotaphicols = c(0,0.15,0.3,0.45,0.6)
iotaphilims = c(0,0.66)
ellomegacols = c(0.6,0.8,1.0,1.2,1.4)
ellomegalims = c(0.6,1.41)

q1 <- ggplot() +
geom_raster(data=iotaVell,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") +
scale_fill_manual(labels=c("limit cycle","1 stable","2 stable","crash","stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=iotaphilims, breaks=iotaphicols, labels=iotaphicols) +
scale_y_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
ylab(expression(paste("Learning rate, ", '\u2113'))) +
xlab(expression(paste("Inflow rate, ", iota))) +
coord_fixed(ratio = 0.75)
ggsave(q1,filename="coop_iotaVell.png", width = 3.5, height = 3.5)

q2 <- ggplot() +
geom_raster(data=iotaVomega,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") +
scale_fill_manual(labels=c("limit cycle","1 stable","2 stable","crash","stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=iotaphilims, breaks=iotaphicols, labels=iotaphicols) +
scale_y_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
ylab(expression(paste("Outflow rate, ", omega))) +
xlab(expression(paste("Inflow rate, ", iota))) +
coord_fixed(ratio = 0.75)
ggsave(q2,filename="coop_iotaVomega.png", width = 3.5, height = 3.5)

q3 <- ggplot() +
geom_raster(data=iotaVphi,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") +
scale_fill_manual(labels=c("limit cycle","1 stable","2 stable","crash","stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=iotaphilims, breaks=iotaphicols, labels=iotaphicols) +
scale_y_continuous(expand=c(0,0), limits=iotaphilims, breaks=iotaphicols, labels=iotaphicols) +
ylab(expression(paste("Resusceptibility rate, ", phi))) +
xlab(expression(paste("Inflow rate, ", iota))) +
coord_fixed(ratio = 1)
ggsave(q3,filename="coop_iotaVphi.png", width = 3.5, height = 3.5)

q4 <- ggplot() +
geom_raster(data=ellVomega,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") +
scale_fill_manual(labels=c("limit cycle","1 stable","2 stable","crash","stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
scale_y_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
ylab(expression(paste("Outflow rate, ", omega))) +
xlab(expression(paste("Learning rate, ", '\u2113'))) +
coord_fixed(ratio = 1)
ggsave(q4,filename="coop_ellVomega.png", width = 3.5, height = 3.5)

q5 <- ggplot() +
geom_raster(data=phiVell,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") +
scale_fill_manual(labels=c("limit cycle","1 stable","2 stable","crash","stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=iotaphilims, breaks=iotaphicols, labels=iotaphicols) +
scale_y_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
ylab(expression(paste("Learning rate, ", '\u2113'))) +
xlab(expression(paste("Resusceptibility rate, ", phi))) +
coord_fixed(ratio = 0.75)
ggsave(q5,filename="coop_phiVell.png", width = 3.5, height = 3.5)

q6 <- ggplot() +
geom_raster(data=phiVomega,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="none") +
scale_fill_manual(labels=c("limit cycle","1 stable","2 stable","crash","stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=iotaphilims, breaks=iotaphicols, labels=iotaphicols) +
scale_y_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
ylab(expression(paste("Outflow rate, ", omega))) +
xlab(expression(paste("Resusceptibility rate, ", phi))) +
coord_fixed(ratio = 0.75)
ggsave(q6,filename="coop_phiVomega.png", width = 3.5, height = 3.5)

## Function to extract legend
g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

legend <- ggplot() +
geom_raster(data=iotaVell,aes(x=V1,y=V2,fill=factor(V3))) +
theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),legend.position="bottom") +
scale_fill_manual(labels=c("0"="cycle","1"="1 stable","2"="2 stable","10"="crash","11"="stable & crash"),values=cols) +
scale_x_continuous(expand=c(0,0), limits=iotaphilims, breaks=iotaphicols, labels=iotaphicols) +
scale_y_continuous(expand=c(0,0), limits=ellomegalims, breaks=ellomegacols, labels=ellomegacols) +
ylab(expression(paste("Learning rate, ", '\u2113'))) +
xlab(expression(paste("Inflow rate, ", iota))) +
coord_fixed(ratio = 0.75) +
guides(fill=guide_legend(title="Nature of equilibria: "))

legend <- g_legend(legend)
ggsave(legend,filename="coop_legend.png", width = 10.5, height = 3.5)
"""
