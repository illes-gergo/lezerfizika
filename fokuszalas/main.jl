using Plots, DelimitedFiles, FFTW, LazyGrids, Base.Threads

include("fuggvenyek.jl")


f = 5e-2;
w0 = 1e-2;
t = 0
z = range(-1e-2, 1e-2, 5000) .+ f
z_ = range(-1e-2, 1e-2, 5000) .+ f
x = range(-1e-2, 1e-2, 1000) .* 2 

Z, X = ndgrid(z_, x)

plotlyjs()
rawData = readdlm(pwd() * "/fokuszalas/ATHz.txt", ',')
nu = rawData[:, 1] * 1e12
spektrum = rawData[:, 2]
#display(plot(spektrum))

@time eredmeny = superPos(2, 1000, nu, spektrum, z, 0)
display(plot(real(eredmeny)))

@time eredmeny2D = superPos(2, 1000, nu, spektrum, Z, X)
#eredmeny2D = dropdims(eredmeny2D, dims=3)
display(heatmap(real(eredmeny2D)))
