using Plots, DelimitedFiles, FFTW, LazyGrids, Base.Threads, DSP

include("fuggvenyek.jl")
default(size=[800, 600])

f = 5e-2;
w0 = 1e-2;
t = 0
mult = 1.0
z = range(-1e-2, 1e-2, 5000) .* 2 .+ mult * f
z_ = range(-1e-2, 1e-2, 5000) .* 2 .+ mult * f
x = range(-1e-2, 1e-2, 1000) .* 2

Z, X = ndgrid(z_, x)

plotlyjs()
rawData = readdlm(pwd() * "/fokuszalas/FocusData.txt", ',')
nu = rawData[:, 1]
spektrum = rawData[:, 2]
fazis = unwrap(rawData[:, 3])

#display(plot(fazis))

 @time eredmeny = superPos(2, 1000, nu, spektrum, z, 0, t, fazis)
display(plot(real(eredmeny)))


@time eredmeny2D = superPos(2, 1000, nu, spektrum, Z, X, t, fazis)
#eredmeny2D = dropdims(eredmeny2D, dims=3)
display(heatmap(real(eredmeny2D)))

