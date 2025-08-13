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

