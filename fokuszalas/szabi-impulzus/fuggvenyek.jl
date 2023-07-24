function gaussian_beam(z, x, t, A, w0, omega, phi, c)

    k_j = k_j_(omega, c)
    z_0_j = z_0_j_(k_j, w0)
    w_j = w_j_(z, w0, z_0_j)

    return A .* w0 ./ w_j .* exp.(-x .^ 2 ./ w_j .^ 2) .* sin.(omega .* t .- k_j .* z .- (k_j .* x .^ 2 .* z) ./ (2 .* (z .^ 2 .+ z_0_j .^ 2)) + atan.(z ./ z_0_j) .+ phi)
end

function gaussian_beam_abs(z, x, t, A, w0, omega, phi, c)

    k_j = k_j_(omega, c)
    z_0_j = z_0_j_(k_j, w0)
    w_j = w_j_(z, w0, z_0_j)

    return A .* w0 ./ w_j .* exp.(-x .^ 2 ./ w_j .^ 2) .* exp.(1im.*(omega .* t .- k_j .* z .- (k_j .* x .^ 2 .* z) ./ (2 .* (z .^ 2 .+ z_0_j .^ 2)) + atan.(z ./ z_0_j) .+ phi))
end

function w_j_(z, w0, z_0_j)
    return w0 .* sqrt.(1 .+ z .^ 2 ./ z_0_j .^ 2)
end

function z_0_j_(k_j, w0)
    return 1 ./ 2 .* k_j .* w0 .^ 2
end

function k_j_(omega_j, c)
    omega_j ./ c
end

function supPos(z, x, inpData, t, wshift, zshift, tshift)
    ret_val = zeros(size(z)..., nthreads())
    central_frq = sum(inpData[:, 2] .* inpData[:, 1]) ./ sum(inpData[:, 2])
    @threads for i in axes(inpData, 1)
        ret_val[:, :, threadid()] .+= gaussian_beam(z .+ (zshift .* (inpData[i, 1] .- central_frq)), x .+ (wshift .* (inpData[i, 1] .- central_frq)), t .+ (tshift .* (inpData[i, 1] .- central_frq)), inpData[i, 2], inpData[i, 3], 2 * pi * inpData[i, 1], inpData[i, 4], 3e8)
    end
    return dropdims(sum(ret_val, dims=3), dims=3)
end

function supPosABS(z, x, inpData, t, wshift, zshift, tshift)
    ret_val = Complex.(zeros(size(z)..., nthreads()))
    central_frq = sum(inpData[:, 2] .* inpData[:, 1]) ./ sum(inpData[:, 2])
    @threads for i in axes(inpData, 1)
        ret_val[:, :, threadid()] .+= gaussian_beam_abs(z .+ (zshift .* (inpData[i, 1] .- central_frq)), x .+ (wshift .* (inpData[i, 1] .- central_frq)), t .+ (tshift .* (inpData[i, 1] .- central_frq)), inpData[i, 2], inpData[i, 3], 2 * pi * inpData[i, 1], inpData[i, 4], 3e8)
    end
    return dropdims(abs.(sum(ret_val, dims=3)), dims=3)
end