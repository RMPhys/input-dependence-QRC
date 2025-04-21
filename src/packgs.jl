# Set the number of BLAS (MKL) threads
""" CAREFUL! changing the number of threads changes 
    the seed. The results of the paper were obtained 
    with ENV["MKL_NUM_THREADS"] = "2"
    while the figure of CPU running time was obtained 
    with ENV["MKL_NUM_THREADS"] = "1"   
"""
ENV["MKL_NUM_THREADS"] = "2"  
using MKL #It speeds up the code a bit
using LinearAlgebra

