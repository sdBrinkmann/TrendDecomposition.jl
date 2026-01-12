println("Testing...")

using Random, Distributions, Statistics, Test
using TrendDecomposition

Random.seed!(7886)

ϵ = rand(Normal(0, 16), 200)

ζ = rand(Normal(0, 1), 200)

μ = zeros(200)
β = zeros(200)


for t in 2:200
    β[t] = β[t-1] + ζ[t]
end

for t in 2:200
    μ[t] = μ[t-1] + β[t-1]
end



x = range(1, 200, 200)

y = μ + ϵ


q = 1 / 16^2
λ = 1 / q


@testset "Iteration hpFilter" begin
    @test hpFilter(y, λ) ≈ hpFilter(y, λ, 1)
    #@test HP(HP(y, λ), λ) == HP(y, λ, 2)

end

@testset "Equivalence hpFilter" begin
    @test hpFilter(y, 0) == y
    @test hpFilter(y, λ) == bohlmannFilter(y, 2, λ)
    @test bhpFilter(y, λ) ≈ hpFilter(y, λ)
end

@testset "Errors hpFilter" begin
    @test_throws AssertionError hpFilter([1], λ)
    @test_throws AssertionError bhpFilter([1], λ)
    @test_throws AssertionError bohlmannFilter([1], 1, λ)
end


## Testing Moving Averages

x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

weights = [.2, .1, .7]


@testset "Equivalence rollingAverage" begin
    @test rollingAverage(x, 3, weights, centered = true, offset=-1) ==
        rollingAverage(x, 3, weights, centered=false)
     @test rollingAverage(x, 3, centered = true, offset=-1) ==
        rollingAverage(x, 3, centered=false)
    @test rollingAverage(x, 5, [.2, .2, .2, .2, .2]) ≈ rollingAverage(x, 5)
end

@testset "Result correctness rollingAverage" begin
    @test rollingAverage(x, 5)[1] == sum(x[1:3]) / 3
    @test rollingAverage(x, 3)[5] == sum(x[4:6]) / 3
    @test rollingAverage(x, 3)[10] == sum(x[9:10]) / 2
    @test rollingAverage(x, 5, centered = false)[10] == sum(x[6:10]) / 5
    @test rollingAverage(x, 5, centered = false, offset=1)[10] == sum(x[7:10]) / 4
    @test rollingAverage(x, 5, offset=-2)[1] == x[1]
    @test rollingAverage(x, 5, offset=2)[5] == sum(x[5:9]) / 5
    
    @test rollingAverage(x, 3, weights)[1] == sum(x[1:2], weights[2:3] / sum(weights[2:3]))
    @test rollingAverage(x, 3, weights)[5] == sum(x[4:6], weights)
end


@testset "Errors rollingAverage" begin
    @test_throws AssertionError rollingAverage(x, 20)
    @test_throws AssertionError rollingAverage(x, 5, offset=5)
    @test_throws AssertionError rollingAverage(x, 5, offset=-5)
    @test_throws AssertionError rollingAverage(x, 3, weights, offset=5)
    @test_throws AssertionError rollingAverage(x, 3, weights, centered=false, offset=-1)
    @test_throws AssertionError rollingAverage(x, 3, centered=false, offset=-1)
    @test_throws AssertionError rollingAverage(x, 4, weights, offset=5)
    @test_throws AssertionError rollingAverage(x, 2, weights, offset=5)
end

@testset "maSeason" begin
    @test maSeason([1, 2], 2) == [1, 2]
    @test maSeason([1, 2, 1, 2, 1, 2, 1, 2], 2) == [1, 2]
    @test maSeason([1, 2, 1, 2, 1, 2, 1, 2], 2, repeating=true) == [1, 2, 1, 2, 1, 2, 1, 2]
    @test maSeason([1, 2, 3, 4, 1, 2, 3, 4], 4) == [1, 2, 3, 4]
end


@testset "EWMA" begin
    @test all(x -> x == y[2] - y[1], holtLinear(y, 0.5, 0.0)[2][2:end, 2])
    @test holtLinear(y, 1.0, 0.5)[2][2:end, 1] == y[2:end]

    @test all(x -> x == y[2] - y[1], holtWinters(y, 0.5, 0.0, 0.5, 1)[2][2:end, 2])
    @test holtWinters(y, 1.0, 0.5, 0.5, 1)[2][2:end, 1] == y[2:end]

    @test brownLinear(y, 1.0; φ = 0.8)[2][3, 2] == (y[2] - y[1]) * 0.8
    @test holtLinear(y, .5, 0.0, φ = 0.8)[2][3, 2] == (y[2] - y[1]) * 0.8  
end

@testset "expSmooth" begin
    sm = expSmoothing([1, 2, 3], 0.5)
    cal = ones(3)
    cal[1] = 0.5
    cal[2] = 0.5 * cal[1] + 0.5 * 2
    cal[3] = 0.5 * cal[2] + 0.5 * 3
    @test cal == sm
end

@testset "EWMA forecasting" begin
    f = holtLinear(y, .84, .84, h = 4) 
    @test f[1][(end-3):end] == [f[2][end, 1] + j * f[2][end, 2] for j in 1:4]

    f = brownLinear(y, .84, h = 4) 
    @test f[1][(end-3):end] == [f[2][end, 1] + j * f[2][end, 2] for j in 1:4]

    f = holtWinters(y, .84, .84, .84, 4, h = 4)
    S = f[2][(end-3):end, 3]
    @test f[1][(end-3):end] == [f[2][end, 1] + j * f[2][end, 2] + S[j] for j in 1:4]

    f = holtWinters(y, .84, .84, .84, 4, h = 4, model=:mul)
    S = f[2][(end-3):end, 3]
    @test f[1][(end-3):end] == [(f[2][end, 1] + j * f[2][end, 2]) * S[j] for j in 1:4]
end



@testset "Holt Brown equivalence" begin
    ω = 0.8

    λ1 = 1 - ω^2
    λ2 = (1 - ω) / (1 + ω)

    res1 = brownLinear(y, ω)
    res2 = holtLinear(y, λ1, λ2)

    @test all(res1[1][3:end] .≈ res2[1][3:end])
    @test all(res1[2][3:end,1] .≈ res2[2][3:end,1])
    @test all(res1[2][3:end,2] .≈ res2[2][3:end,2])

    res3 = brownLinear(x, ω, startValues = (.5, .5))
    res4 = holtLinear(x, λ1, λ2, startValues = (.5, .5))

    @test all(res3 .≈ res4)   
end


@testset "Optimization" begin
    f, D, λ = holtWinters(y, 4)
    f1, D1 = holtWinters(y, λ[1], λ[2], λ[3], 4, φ = λ[4])
    @test f[3:end] == f1[3:end]
    @test D[3:end, :] == D1[3:end, :]

    f, D, λ = holtWinters(y, 4, φ = 1.0)
    f1, D1 = holtWinters(y, λ[1], λ[2], λ[3], 4)
    @test f[3:end] == f1[3:end]
    @test D[3:end, :] == D1[3:end, :]

    f, D, λ = holtLinear(y)
    f1, D1 = holtLinear(y, λ[1], λ[2], φ = λ[3])
    @test f[3:end] == f1[3:end]
    @test D[3:end, :] == D1[3:end, :]

    f, D, λ = holtLinear(y, φ = 1.0)
    f1, D1 = holtLinear(y, λ[1], λ[2])
    @test f[3:end] == f1[3:end]
    @test D[3:end, :] == D1[3:end, :]

    f, D, λ = brownLinear(y)
    f1, D1 = brownLinear(y, λ[1], φ = λ[2])
    @test f[3:end] == f1[3:end]
    @test D[3:end, :] == D1[3:end, :]

    f, D, λ = brownLinear(y, φ = 1.0)
    f1, D1 = brownLinear(y, λ)
    @test f[3:end] == f1[3:end]
    @test D[3:end, :] == D1[3:end, :]

end
