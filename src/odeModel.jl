function gain(u,r,K,ind)
    return r[ind] .* u[ind] .* (K[ind] .- u[ind])
end

function holling2(u,r,kappa,A,ind)
    return r[ind] .* u[ind] .* (A' * u)[ind] ./ (kappa[ind] .+ (A' * u)[ind])
end

# function dependentLoss(u,r,kappa,A,ind)
#     L = zeros(size(ind,1))
#     for i = ind
#         L[i - ind[1] + 1] = sum((A[i,:] .* u[i] .* u) ./ (kappa .+ sum(A .* u, dims = 1)[:]))
#     end
#     return L
# end

function dependentLoss(u,r,kappa,A,ind)
    return u[ind] .* sum((A[ind,:]' .* u) ./ (kappa .+ sum(A .* u, dims = 1)[:]), dims = 1)[:]
end

function webDynamic!(du,u,p,t)
    r,K,A,kappa,m,ind = p
    du[ind[1]] = gain(u,r,K,ind[1]) .- m[ind[1]] .* dependentLoss(u,r,kappa,A,ind[1])
    du[ind[2]] = holling2(u,r,kappa,A,ind[2]) .- m[ind[2]] .* dependentLoss(u,r,kappa,A,ind[2]) .- 0.01 .* u[ind[2]]
    du[ind[3]] = holling2(u,r,kappa,A,ind[3]) .- m[ind[3]] .* u[ind[3]]
end
