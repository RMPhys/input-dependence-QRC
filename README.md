# input-dependence-QRC
The Mathematica notebook for some of the examples, generating Fig.1 and the Julia code to generate Fig. 2 of our paper Input-dependence in quantum reservoir computing.

# How to run the Mathematica notebook
If you don't have a Mathematica license, you can at least use [Wolfram Player](https://www.wolfram.com/player/) to see the notebook and run it. 

# Julia code
The notebook Fig2.ipynb generates Fig.2 of the manuscript. The script `sqmodule.jl` contains the necessary code to perform single qubit reservoir calculations. If new single qubit models are to be tested, it is only necessary to create a new `struct` with the abstract type `SQModel` and the corresponding `dynamics` function for a single time step, adding the new model type to the function `build_SQModel`.