function gaussNyalab(x, z, t, A, lambda0, fazis)

    omega0 = 3e8 / lambda0 * 2 * pi
    k = 2 * pi ./ lambda0
    z0 = pi * w0^2 ./ lambda0
    # zM = f ./ (1 + (f ./ z0) .^ 2)
    z_ = z
    recR = @. (-1 / f + z_ .* (1 / f^2 + 1 ./ z0 .^ 2)) ./ ((1 - z_ ./ f) .^ 2 + (z_ ./ z0) .^ 2)
    recWz = @. sqrt((1 ./ z0) ./ ((1 - z_ ./ f) .^ 2 + (z_ ./ z0) .^ 2) * pi ./ lambda0)
    gouy = @. atan(z_ ./ z0)
    # println(omega0)
    V = @. A * w0 .* recWz .* exp(-x .^ 2 .* recWz .^ 2) .* exp(-1im .* (k .* (z_) - omega0 * t + k .* x .^ 2 ./ 2 .* recR - gouy + fazis))
    return V
end

function superPos(startI, endI, nu, spektrum, zGrid, xGrid::LazyGrids.GridSL, t, fazis)
    supPos = zeros(Complex{Float64}, size(zGrid)..., nthreads())
    @threads for ind in startI:endI
        supPos[:, :, threadid()] += gaussNyalab(xGrid, zGrid, t, spektrum[ind], 3e8 / nu[ind], fazis[ind])
    end
    return dropdims(sum(supPos, dims=3), dims=3)
end


function superPos(startI, endI, nu, spektrum, zGrid, xGrid::Number, t, fazis)
    supPos = zeros(Complex{Float64}, size(zGrid))
    for ind in startI:endI
        supPos += gaussNyalab(xGrid, zGrid, t, spektrum[ind], 3e8 / nu[ind], fazis[ind])
    end
    return supPos
end