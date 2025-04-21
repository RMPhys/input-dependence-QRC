"""
    Common

Author: Rodrigo Martínez-Peña, rodrigo.martinez@dipc.org 
This module contains all the functions that are common to the QRC models. 
"""

module Common
include("packgs.jl")
using Statistics # For the cor function.

abstract type Metric end

@kwdef struct Metric_C <: Metric
    name::String = "C"
end 

function metric_val(metric::Metric_C,target,prediction)
    return (cor(target,prediction)[1,1])^2
end


end #end of module


