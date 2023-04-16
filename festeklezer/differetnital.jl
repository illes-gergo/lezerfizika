function diffegy(t, current_values)
    # Az értékek sorrendje n(t), q(t), na(t), tauc(t)
    n = current_values[1]
    q = current_values[2]
    na = current_values[3]
    tauc = max(eta * L^3 / 8 / c / pi^2 * (n * (sigmaell - sigmalll) * V)^2, eta * L / 10 / c)

    dndt = I(t) * sigmap * (N - n) - sigmaell * c / eta * n * q - n / tau

    dqdt = (sigmaell - sigmalll) * c / eta * n * q - q / tauc + omega * n / tau - sigmaaal * c / eta * na * q +
           sigmaeal * c / eta * (Na - na) * q - sigmalal * c / eta * (Na - na) * q

    dnadt = -I(t) * sigmapa * na + sigmaeal * c / eta * (Na - na) * q + (Na - na) / taua - sigmaaal * c / eta * na * q +
            2 * omegaa / taua / L / (sigmaeaa - sigmalaa) * (exp((sigmaeaa - sigmalaa) * (Na - na) * L) - 1)

    return [dndt, dqdt, dnadt]
end

function tauc_(n)
    return max(eta * L^3 / 8 / c / pi^2 * (n * (sigmaell - sigmalll) * V)^2, eta * L / 10 / c)
end

function RK4(f::Function, t, y, step)
    k1 = f(t, y)
    k2 = f(t + step / 2, y + step / 2 * k1)
    k3 = f(t + step / 2, y + step / 2 * k2)
    k4 = f(t + step, y + step * k3)
    return y .+ step .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4) ./ 6
end