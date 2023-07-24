using Plots, DelimitedFiles, LazyGrids, Base.Threads

gr()

c = 3e8
x = range(-10,10,length=200).*1e-3
z = copy(x)

wshift = 10e-15 
zshift = -1e-15
tshift = 0e-24
t = 1.3e-3/c

Z,X = ndgrid(z,x)

include("fuggvenyek.jl")

inpData = readdlm("g0_m1_300GHz_1mm_s3.85_521MVm_300.txt")

central_frq = sum(inpData[:, 2] .* inpData[:, 1]) ./ sum(inpData[:, 2])

#@time supData=supPos(Z,X,inpData,t,wshift,zshift,tshift)

@time supDataABS=supPosABS(Z,X,inpData,t,wshift,zshift,tshift)

#display(heatmap(supData,colormap=:jet))
display(heatmap(supDataABS,colormap=:jet))