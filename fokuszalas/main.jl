using Plots, DelimitedFiles, FFTW

include("fuggvenyek.jl")


f = 5e-2 / 2;
w0 = 1e-2;
t = 0
z = range(-1e-2,1e-2,5000) .+ f

plotlyjs()
rawData = readdlm(pwd() * "/fokuszalas/ATHz.txt", ',')
nu = rawData[:, 1] * 1e12
spektrum = rawData[:, 2]
#plot(nu, spektrum)

@time eredmeny = superPos(2,10000,nu,spektrum,z,0)
plot(real(eredmeny))