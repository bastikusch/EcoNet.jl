## types

struct FoodWeb
    links::Matrix
    nicheValues::Vector{Float64}
    distribution::Vector{Int64}
    index::Vector{UnitRange}
end

## functions

function generateFoodWeb(S::Int, C::Float64)
    @assert (C <= 0.5 && C >= 0.0)
    β = 1.0 / (2.0 * C) - 1.0
    n = vcat(0.0,sort(rand(S - 1)))
    # n = sort(rand(S))
    r = n .* rand(Beta(1,β),S)
    c = (n .- (r ./ 2)) .* rand(S) .+ (r ./ 2)
    A = getAdjacency(n, r, c)
    while (isolatedSpecies(A))
        replaceIsolatedSpecies!(A, n, r, c, β)
        A = getAdjacency(n, r, c)
    end
    T = findall(sum(A, dims=2)[:] .== 0)
    B = findall(sum(A, dims=1)[:] .== 0)
    I = findall((sum(A, dims=2)[:] .> 0) .& (sum(A, dims=1)[:] .> 0))
    newInd = vcat(B,I,T)
    distri = [size(B,1), size(I,1), size(T,1)]
    ind = [1:distri[1],distri[1]+1:size(B,1) + size(I,1),size(B,1) + size(I,1) + 1:S]
    A = getAdjacency(n[newInd], r[newInd], c[newInd])
    return FoodWeb(A,n[newInd],distri,ind)
end

function isolatedSpecies(A::Matrix)
    return size(findall((sum(A, dims=2)[:] .== 0) .& (sum(A, dims=1)[:] .== 0)),1) > 0
end

function replaceIsolatedSpecies!(A, n, r, c, β)
    ind = findall((sum(A, dims=2)[:] .== 0) .& (sum(A, dims=1)[:] .== 0))
    # ind = findall(x -> x != 1, ind)
    for i in ind
        n[i] = rand()
        r[i] = n[i] * rand(Beta(1,β))
        c[i] = (n[i] - (r[i] / 2)) * rand().+ (r[i] / 2)
    end
end

function getAdjacency(n::Vector,r::Vector,c::Vector)
    A = zeros(size(n,1),size(n,1))
    for i = 1:size(n,1)
        for j = 2:size(n,1)
            A[i,j] = (n[i] <= c[j] + r[j] / 2) && (n[i] >= c[j] - r[j] / 2)
        end
    end
    return A - Diagonal(A)
end
