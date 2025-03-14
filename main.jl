#import Pkg
#Pkg.add("DataFrames")
#Pkg.add("CSV")
#Pkg.add("Plots")
#Pkg.add("SpecialFunctions")
#Pkg.add("QuadGK")
using DataFrames
using CSV
using Plots
using SpecialFunctions
using QuadGK

zer = CSV.File("zeroes.csv")
#to run in vscode uncomment the import pkg and the pkg add whatever and then right click and select run code
function usefulEvalPoints(startVal, endVal)
    output = round.(startVal):round(endVal)
    return vcat(output .+ 0.01, output .- 0.01)
end
q = 0.05
startVal = 2.0
endVal = 10.0
t = sort(vcat(range(startVal, endVal, step=q), usefulEvalPoints(startVal, endVal)))
c = 100

function T(x)
    return Li(x) .- sumg(x)
end

function sumg(x)
    output2 = fill(0.0, length(t))
    for n in 1:c
        output2 .+= g(x, zeroes(n))
    end
    return 2 .* output2
end

function g(x, b)
    output3 = fill(0.0, length(t))
    for k in 1:20
        #b, k scalar x arr
        output3 .+= ((mobius(k) / k) .* integ(x, b, k))
    end
    return output3
end

function zeroes(n)
    return CSV.getcolumn(zer, 1).column[n]
end

function integ(tVals, b, n)
    #f(d) = map(z -> exp(z + (b * im * l(x)))/(z + (b * im * l(x))), d)
    f_b(x) = (b * 1im)*log(x)/n
    f(x) = real(quadgk(z -> exp(z + f_b(x))/(z + f_b(x)), -Inf, 0.5*log(x)/n)[1])
    return map(x -> f(x), tVals)
end

function Li(x)
    f(d) = map(p -> 1/log(p), d)
    return map(i -> quadgk(f, 2, i)[1], x)
end

function mobius(n)
    arr = coprime(n)
    arr = exp.(2im * π * arr / n)
    return round(real(sum(arr)))
end

function coprime(n)
    output = 1:n
    return output[gcd.(output, n) .== 1]
end

function startProgram()
    start = time_ns()
    s = T(t)
    println(time_ns() - start)
    #anim = @animate for i in 1:100
    #    plot(t, T f(t, i))
    #end
    #gif(anim, fps=15)
    #@gif for i in 1:20
    display(plot(t, s, xticks = 0:5:100))
    #end
    readline()
    #grid(true)
    #display(plot!)
end

startProgram()
