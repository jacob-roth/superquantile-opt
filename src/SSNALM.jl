using LinearAlgebra
using Random
using Printf
using SparseArrays
using StatsBase
using Dates, DataStructures, DelimitedFiles, LaTeXStrings, LaTeXTabulars # save results
using Gurobi, SparseMatricesCSR, SparseArrays, OSQP # experiments
if Sys.iswindows()
  global const PROJPATH = match(r".*cc-non", @__DIR__).match*"\\"
else
  global const PROJPATH = match(r".*cc-non/", @__DIR__).match
end

include("structures.jl")
include("projection.jl")
include("obj.jl")
include("genprob.jl")
include("utility.jl")
include("diagnostics.jl")
include("phi.jl")
include("ssn.jl")
include("alm.jl")
include("linesearch.jl")
include("proj.jl")
include("admm.jl")
include("cvar_solver.jl")
include("cvar_oa_solver.jl")
include("qr_solver.jl")
