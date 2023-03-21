using Base.Threads, Plots

include("differetnital.jl")

N = 9e24 # 1/m^3
Na = 1.08e23 # 1/m^3
tau = 4e-9
taua = 4e-9
sigmap = 1.15e-21
sigmapa = 3.8e-21

sigmaell = 1.1e-20
sigmalll = 0.4e-20
sigmaaal = 4.1e-20
sigmaeal = 0.9e-20
sigmalal = 1.3e-20
sigmaeaa = 3.05e-20
sigmalaa = 2e-20

c = 3e8
eta = 1.32
omega = 5.2e-8
omegaa = 2.5e-6
L = 0.3e-2
V = 0.38

Imax = 1e20
FWHM = 9e-9

I(t) = Imax * exp(-t^2 / FWHM^2)

dt = 1e-13

tspan = 20e-9

t = range(-tspan / 2, tspan / 2, step=dt)

calculatedValues = Array{Float64}(undef, 3, length(t))

calculatedValues[1, 1] = N
calculatedValues[2, 1] = 0
calculatedValues[3, 1] = Na

for ii in 1:(length(t)-1)
    calculatedValues[:, ii+1] = RK4(diffegy, t[ii+1], calculatedValues[:, ii], dt)
end



display(plot(calculatedValues[1,:]))
display(plot(calculatedValues[2,:]))
display(plot(calculatedValues[3,:]))