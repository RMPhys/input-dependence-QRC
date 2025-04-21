"""
    SQmodule

Author: Rodrigo Martínez-Peña, rodrigo.martinez@dipc.org 
This module contains all the functions about the Single Qubit method that
we will employ in the main. 
"""
module SQmodule
include("packgs.jl")
using SparseArrays
using LinearAlgebra
const σz = sparse(Matrix{ComplexF64}([1 0; 0 -1])) 
const σx = sparse(Matrix{ComplexF64}([0 1; 1 0]))
const σy = sparse(Matrix{ComplexF64}([0 -1im; 1im 0]))

abstract type SQModel end

function Ry(s::Float64)
    return Matrix{ComplexF64}([cos(s/2.0) -sin(s/2.0); sin(s/2.0) cos(s/2.0)]) 
end

function Rx(s::Float64)
    return Matrix{ComplexF64}([cos(s/2.0) -1im*sin(s/2.0); -1im*sin(s/2.0) cos(s/2.0)]) 
end

function Rz(s::Float64)
    return Matrix{ComplexF64}([exp(-1im*s/2.0) 0; 0 exp(1im*s/2.0)])
end

function build_SQModel(tupla::Tuple{Vararg{Float64}},type::Symbol)
    if type == :Example24
        return Example24(tupla...)
    else
        error("type must be :Example24")
    end
end

#########################Example24################################

struct Example24 <: SQModel 
    g :: Float64 #input strength
    ϵ :: Float64 #reset rate
    λ :: Float64 #regularization parameter
    σ :: Matrix{ComplexF64} #reset matrix
    name :: String #name of the model
    function Example24(g :: Float64,ϵ :: Float64,λ :: Float64)
        σ = Matrix{ComplexF64}([1 0; 0 0])
        return new(g,ϵ,λ,σ,"Example24")
    end
end   
function dynamics(model::Example24,ρin::Matrix{ComplexF64},uin::Vector{Float64})
    U = Ry(model.g*uin[1])
    ρprima = U * ρin * U' # encoding channel
    ρout = (1.0-model.ϵ) * ρprima + model.ϵ * model.σ # contracted channel
    return ρout
end 

####################################################################


function washout_phase(model::SQModel,u::Matrix{Float64},washout::Int64)
    ρ0 = zeros(ComplexF64,2,2)
    ρ0[1,1] = 1.0+0.0im
    # Washing out phase
    for t=1:washout
        ρ0  = dynamics(model,ρ0,u[t,:])
    end
    return ρ0
end


function RES_construct(model::SQModel,u::Matrix{Float64},ρ0::Matrix{ComplexF64},N::Int64)
    RES = Matrix{Float64}(undef,N,3)
    # Training phase
    ρ = 1*ρ0
    for t=1:N
       ρ = dynamics(model,ρ,u[t,:])
       RES[t,1] = real(tr(σx*ρ))
       RES[t,2] = real(tr(σy*ρ))
       RES[t,3] = real(tr(σz*ρ))
    end 
    return RES,ρ
end

function train(model::SQModel,u::Matrix{Float64},y::Matrix{Float64},N::Int64,washout::Int64)
    # Washing out phase
    ρ0 = washout_phase(model,u[1:washout,:],washout)
    # Training phase
    RES,ρ = RES_construct(model,u[washout+1:N,:],ρ0,N-washout)
    RESplusone = hcat(RES,ones(Float64,size(RES, 1)))
    weights = regression(RESplusone,y[washout+1:N,:],model.λ) 
    return weights,ρ
end


function pred_closedloop(model::SQModel,u1::Vector{Float64},weights::Matrix{Float64},ρ0::Matrix{ComplexF64},Ntest::Int64)   
    l = size(weights,2)
    yhat = Matrix{Float64}(undef,Ntest,l)
    x = Vector{Float64}(undef,3)
    ρ = 1*ρ0 
    u = 1*u1
    # Next steps
    for t=1:Ntest
        ρ = dynamics(model,ρ,u)
        x[1] = real(tr(σx*ρ))
        x[2] = real(tr(σy*ρ))
        x[3] = real(tr(σz*ρ))
        yhat[t,:] = transpose(weights[1:end-1,:])*x .+ weights[end,:]
        u = yhat[t,:]
    end
    return yhat
end



function pred_openloop(model::SQModel,u::Matrix{Float64},weights::Matrix{Float64},ρ0::Matrix{ComplexF64},Ntest::Int64) 
    RES,ρ = RES_construct(model,u,ρ0,Ntest)
    RESplusone = hcat(RES,ones(Float64,size(RES, 1)))
    yhat=RESplusone*weights  
    return yhat 
end


function regression(A::Matrix{Float64},y::Matrix{Float64},λ::Float64) 
    weights = (transpose(A)*A + λ*I)\(transpose(A)*y)
    return weights
end



end #end of module
