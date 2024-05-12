# using Profile, Serialization
import Pkg; Pkg.activate("."); Pkg.instantiate()
# using MKL
using JuMP
include("../../src/SSNALM.jl")
include("experiments.jl")
# GC.enable_logging(true)
BLAS.get_num_threads()
# BLAS.set_num_threads(1)

#
# initial parameters
#

# parameters
Tf = Float64
Ti = Int64
Vf = AbstractVector{Tf}
Vi = AbstractVector{Ti}

# alm parameters
lin_solver = (:M, 1)
AP = ALMParameters()
AP.maxait = 100
AP.stoptol = 1e-6

# ssn parameters
SP = SSNParameters(lin_solver)
SP.maxnit = 200
SP.maxlit = 50
SP.tolsub = 1e-7
SP.ls_tau = 0.7
SP.ls_c1 = 1e-4
SP.ls_c2 = 0.9
SP.ls_method = :sw # {:bt, :sw, :ww}

# options
options = default_options()

#
# m scaling and linear solver
#

# include("experiment_m_linsolver.jl")

#
# n & m scaling
#

DATAPATH = "/home/roth0674/drive/nm_5pct"
include("experiment_nm.jl")

#
# diagnose
#

out_n, TI_n = reconstruct_output(DATAPATH);
TI = TI_n[11,2,3]
# TI = TI_n[3,2,3]
# TI = TI_n[3,2,3]

CI = ConstraintInstance(TI.seed, TI.n, TI.m, get_k(TI), TI.scale_row, TI.scale_col, TI.t);
# sum_{i=1:k}(rl[i]) easy: [24843.331704251614, 24615.259224878806]
# sum_{i=1:k}(rl[i]) hard: [22832.393208322253, 22728.26722848328]
# sum_{i=1:k}(rl[i]) used: [22832.393208322253, 22728.26722848328]
OI = ObjectiveInstance(TI.seed, TI.objtype, TI.n);
PI = ProblemInstance(CI, OI);

obj = FunctionValue(Tf, TI.n, TI.objtype); # 2==diagonal
obj.H .= PI.C;

TI.SP.ls_tau = 0.7
TI.SP.ls_c1 = 1e-4
TI.SP.ls_c2 = 0.9
# TI.SP.ls_method = :bt

x0 = zeros(Tf, TI.n);
lambda0 = [zeros(TI.m[l]) for l in eachindex(TI.m)];
mux0 = zeros(Tf, TI.n);
debug = false;
TI.AP.maxait = 220#20 # 4=just before linear subproblem
TI.SP.maxnit = 200
TI.SP.lin_solver = (:n,1)
res0 = alm(x0, lambda0, mux0, TI.sigma0, TI.tau0, PI, TI.AP, TI.SP, debug)
x0, lambda0, mux0 = res0[:x], res0[:lambda], res0[:mu]
TI.SP.maxnit = 100000
TI.AP.maxait = 100
res = alm(x0, lambda0, mux0, TI.sigma0, TI.tau0, PI, TI.AP, TI.SP, debug)

sum(res0[:ctime_nd]) + sum(res0[:ctime_proj]) + sum(res0[:ctime_Atlam]) + sum(res0[:ctime_sort]) + sum(res0[:ctime_diagnostics]) + sum(res0[:ctime_H]) + sum(res0[:ctime_g]) + sum(res0[:ctime_v]) + sum(res0[:ctime_Ax])

nbuff = 7
ncycl = 4
buff = zeros(nbuff)
# buffvals = [rand(10)./(collect(1:10).^2); 0.123; 0.456; 0.789; 0.1231; 0.4561; 0.7891; 0.12311; 0.45611; 0.78911; 0.123111]
buffvals'
for i in 1:length(buffvals)
  idx = mod1.(i-ncycl:i-1,nbuff)
  buff[mod1(i,nbuff)] = buffvals[i]
  if i > nbuff
    println("---------- i=$i ----------")
    println(buff)
    println(buffvals[i-nbuff+1:i])
    @assert(norm(sum(buff)-sum(buffvals[i-nbuff+1:i])) < 5eps())
    num = buff[mod1(i,nbuff)]
    den = buff[mod1.(i-ncycl:i-1,nbuff)]
    v = abs.(1 .- (num ./ den))
    println(idx)
    println("num = $num")
    println("den = $den")
    println("v  = $v\n")
  end
end

buffvals = [
+9.863e+00 
+9.731e+00 
+7.890e+00 
+8.984e+00 
+9.009e+00 
+8.824e+00 
+8.200e+00 
+8.005e+00 
+7.517e+00 
+7.257e+00 
+2.343e+00 
+3.243e+00 
+3.519e+00 
+3.692e+00 
+3.681e+00 
+3.591e+00 
+3.315e+00 
+3.102e+00 
+2.699e+00 
+2.551e+00 
+2.508e+00
]


resnrm=[
 +2.395e-03
 +1.223e-03
 +6.114e-04
 +6.608e-04
 +2.395e-03
 +1.223e-03
 +6.114e-04
 +6.608e-04
 +2.395e-03
 +1.223e-03
 +6.114e-04
 +6.608e-04
 +2.395e-03
 +1.223e-03
 +6.114e-04
 +6.608e-04
 +2.395e-03
]


#
# test ls params
#

# 1356 lsit, 228 nit, 31sec
TI.SP.ls_c1 = 1e-2 # avoids cycling on one case...
TI.SP.ls_c2 = 2e-2

# 931 lsit, 255 nit, 24sec
TI.SP.ls_c1 = 1e-1 # avoids cycling on one case...
TI.SP.ls_c2 = 2e-1

# 908 lsit, 259 nit, 24sec
TI.SP.ls_c1 = 0.2 # avoids cycling on one case...
TI.SP.ls_c2 = 0.25

# 1097 lsit, 271 nit, 28sec
TI.SP.ls_c1 = 0.4 # avoids cycling on one case...
TI.SP.ls_c2 = 0.45

# 931 lsit, 255 nit, 24sec
TI.SP.ls_c1 = 1e-1 # avoids cycling on one case...
TI.SP.ls_c2 = 0.2
sum(res0[:lsit])
sum(res0[:nit])

# 1362 lsit, 458 nit, 36sec
TI.SP.ls_c1 = 1e-1 # avoids cycling on one case...
TI.SP.ls_c2 = 0.5
sum(res0[:nit])
sum(res0[:lsit])

# 1291 lsit, 384 nit, 34sec
TI.SP.ls_c1 = 1e-1 # avoids cycling on one case...
TI.SP.ls_c2 = 0.9


#? IS THERE A TRADEOFF BETWEEN LSIT AND NIT BASED ON C1 AND C2?
TI.SP.ls_method = :sw # {:bt, :sw, :ww}
c1c2 = [
  [
    #1
    [0.001,0.0011],
    [0.001,0.0021],
    [0.001,0.011],
    [0.001,0.1], # 28 sec
    [0.001,0.9],
  ],
  [
    #2
    [0.01,0.011],
    [0.01,0.021],
    [0.01,0.11],
    [0.01,0.51],
    [0.01,0.91],
  ],
  [
    #3
    [0.1,0.11],
    [0.1,0.21], # 24sec
    [0.1,0.31],
    [0.1,0.51],
    [0.1,0.91],
  ],
  [
    #4
    [0.5,0.51],
    [0.5,0.61],
    [0.5,0.81],
    [0.5,0.91],
    [0.5,0.99],
  ],
  # [
  #   #5: BAD
  #   [0.9,0.901],
  #   [0.9,0.91],
  #   [0.9,0.99],
  #   [0.9,0.999],
  #   [0.9,0.9999]
  # ]
]
# out = Dict()
for exper in eachindex(c1c2)
  out[exper] = zeros(4,length(c1c2[exper])) # ait, nit, lsit, time
  for j in eachindex(c1c2[exper])
    TI.SP.ls_c1 = c1c2[exper][j][1]
    TI.SP.ls_c2 = c1c2[exper][j][2]
    x0 = zeros(Tf, TI.n);
    lambda0 = [zeros(TI.m[l]) for l in eachindex(TI.m)];
    mux0 = zeros(Tf, TI.n);
    debug = false;
    TI.AP.maxait = 120#20 # 4=just before linear subproblem
    TI.SP.maxnit = 200
    TI.SP.lin_solver = (:n,1)
    res0 = alm(x0, lambda0, mux0, TI.sigma0, TI.tau0, PI, TI.AP, TI.SP, debug)
    out[exper][1,j] = sum(res0[:iter])
    out[exper][2,j] = sum(res0[:nit])
    out[exper][3,j] = sum(res0[:lsit])
    out[exper][4,j] = sum(res0[:walltime])
  end
end
ax.cla()
for exper in eachindex(c1c2)
  if exper <= 4
    ax.plot(out[exper][4,:],label=exper)
  end
end
ax.legend()

TI.SP.ls_c1 = 1e-4
TI.SP.ls_c2 = 0.9

#
# integrate proj
#

y0sort = vec(readdlm("y.csv"))
k = Int(readdlm("k.csv")[1])
r = readdlm("r.csv")[1]

# y0sort = Rational{BigInt}.(vec(readdlm("y.csv")))
# k = Int(readdlm("k.csv")[1])
# r = Rational{BigInt}(readdlm("r.csv")[1])

ybarsort1 = similar(y0sort)
ybarsort2 = similar(y0sort)
ybarsort3 = similar(y0sort)
active = sum(y0sort[1:k]) > r


prepop = false
verb = true
spl1, (k0l1, k1l1), _, _, _, lam1 = @views project_maxksum_grid!( ybarsort1, y0sort, r, k, active, prepop, verb)
spl2, (k0l2, k1l2), _, _, _, lam2 = @views project_maxksum_snake!(ybarsort2, y0sort, r, k, active, prepop, verb)
spl3, (k0l3, k1l3), _, _, _, lam3 = @views project_maxksum_lcp!(  ybarsort3, y0sort, r, k, active, prepop, verb)

θ = 0.4629401236789819
λ = 0.049004338058988184

#=
#! comment out
#
# diagnose
#

out_n, TI_n = reconstruct_output(PROJPATH * "results/tune/m_n1solver");
TI = TI_n[11,1,1]
out_n[11,1,1]

TI = TIs[13,1,1]
out[13,1,1]


CI = ConstraintInstance(TI.seed, TI.n, TI.m, get_k(TI), TI.scale_row, TI.scale_col, TI.t);
OI = ObjectiveInstance(TI.seed, TI.objtype, TI.n);
PI = ProblemInstance(CI, OI);

CI.xlo .= -Inf
CI.xhi .= +Inf
CI.xloidx .= false
CI.xhiidx .= false
CI.Xunconstr = true
OI = ObjectiveInstance(TI.seed, TI.objtype, TI.n);
PI = ProblemInstance(CI, OI);

obj = FunctionValue(Tf, TI.n, TI.objtype); # 2==diagonal
obj.H .= PI.C;

x0 = zeros(Tf, TI.n);
lambda0 = [zeros(TI.m[l]) for l in eachindex(TI.m)];
mux0 = zeros(Tf, TI.n);
debug = false;
TI.AP.maxait = 20#20 # 4=just before linear subproblem
TI.SP.maxnit = 100
TI.SP.lin_solver = (:n,1)
res0 = alm(x0, lambda0, mux0, TI.sigma0, TI.tau0, PI, TI.AP, TI.SP, debug)
x0, lambda0, mux0 = res0[:x], res0[:lambda], res0[:mu]
TI.SP.maxnit = 100000
TI.AP.maxait = 100
res = alm(x0, lambda0, mux0, TI.sigma0, TI.tau0, PI, TI.AP, TI.SP, debug)

res = res0
function compare_gurobi_gj(res::Dict, PI::ProblemInstance)
  #=
  define R := A' * (I-G_B) * A = (phi.H - obj.H - sigma(I - G_X)) / sigma
  define h := A * d
  define hbar = G_B(h) = Π'_B(ytilde; h)
  define R_d_g := gurobi computation of R * d = A' * h - A' * hbar
  define R_d_j := julia computation of R * d = A' * h - A' * hbar
  define R_d := R * d
  compare: R_d and R_d_g
  =#
  x = res[:x];
  sigma = res[:sigma];
  lambda = res[:lambda];
  mu = res[:mu];

  utilde = x .- 1/sigma .* mu;
  ubar = similar(x);
  uhat = similar(x);

  ytilde = [PI.A[l] * x + PI.b[l] - 1/sigma .* lambda[l] for l in 1:PI.L];
  ybar_g = [similar(ytilde[l]) for l in 1:PI.L];
  ybar_j = [similar(ytilde[l]) for l in 1:PI.L];
  yhat = [similar(ytilde[l]) for l in 1:PI.L];
  ytmp = [similar(ytilde[l]) for l in 1:PI.L];
  mubarsort_g = [similar(ytilde[l]) for l in 1:PI.L];
  mubarsort_j = [similar(ytilde[l]) for l in 1:PI.L];
  sig_g = [zeros(Int64, length(ytilde[l])) for l in 1:PI.L];
  sig_j = [zeros(Int64, length(ytilde[l])) for l in 1:PI.L];
  siginv_g = [zeros(Int64, length(ytilde[l])) for l in 1:PI.L];
  siginv_j = [zeros(Int64, length(ytilde[l])) for l in 1:PI.L];
  k0bar_g = zeros(Int64, PI.L);
  k1bar_g = zeros(Int64, PI.L);
  k0bar_j = zeros(Int64, PI.L);
  k1bar_j = zeros(Int64, PI.L);
  for l in 1:PI.L
    outl = project_maxksum_gurobi(ytilde[l], PI.r[l], PI.k[l])
    ybar_g[l] .= outl[:xbar]
    mubarsort_g[l] .= outl[:mu]
    
    sig_g[l] .= sortperm(ybar_g[l], rev=true)
    k0bar_l, k1bar_l = get_k0k1(ybar_g[l][sig_g[l]], PI.k[l])
    k0bar_g[l] = k0bar_l
    k1bar_g[l] = k1bar_l

    sig_j[l] .= sortperm(ytilde[l], rev=true)
  end
  active = Bool[1, 0]
  project_B!(ytilde, sig_j,
    ybar_j, yhat, ytmp, # ProjBWorkspace
    k0bar_j, k1bar_j, # diagnostics
    PI.r, PI.k, active # parameters
  )
  for l in 1:PI.L
    mubarsort_j[l][sig_j[l]] .= get_mu(ybar_j[l][sig_j[l]], ytilde[l][sig_j[l]], PI.r[l], PI.k[l])
  end
  for l in 1:PI.L
    siginv_g[l] .= invperm(sig_g[l])
    siginv_j[l] .= invperm(sig_j[l])
  end
  @assert( maximum(norm.(ybar_g-ybar_j, Inf)) <= 1e-4 )
  @assert( maximum(norm.(mubarsort_g-mubarsort_j, Inf)) <= 1e-3 )
  hcat(mubarsort_g[1][sig_g[1]], mubarsort_j[1][sig_j[1]])
  dddd = mubarsort_g[1][sig_g[1]] .- mubarsort_j[1][sig_j[1]]
  dddd[abs.(dddd) .> 1e-4]

  # test differentiability
  fdiff_j = Vector{Dict}(undef, PI.L)# zeros(Bool, PI.L)
  fdiff_g = Vector{Dict}(undef, PI.L)# zeros(Bool, PI.L)
  for l in 1:PI.L
    fdiff_j[l] = test_fdiff(ytilde[l], ybar_j[l], k0bar_j[l], k1bar_j[l], PI.k[l], PI.r[l])
    fdiff_g[l] = test_fdiff(ytilde[l], ybar_g[l], k0bar_g[l], k1bar_g[l], PI.k[l], PI.r[l])
  end
 
  d = randn(PI.n)
  h = [PI.A[l] * d for l in 1:PI.L];
  hbar_g = [zeros(PI.m[l]) for l in 1:PI.L]; # P(Ad)
  hbar_j = [zeros(PI.m[l]) for l in 1:PI.L]; # P(Ad)
  critcone_g = Vector{Dict}(undef, PI.L)
  critcone_j = Vector{Dict}(undef, PI.L)
  for l in 1:PI.L
    outl_j = project_critcone_gurobi(ytilde[l], ybar_j[l], mubarsort_j[l], PI.k[l], k0bar_j[l], k1bar_j[l], PI.r[l], h[l], fdiff_j[l])
    critcone_j[l] = outl_j
    hbar_j[l] .= outl_j[:hbar];
    outl_g = project_critcone_gurobi(ytilde[l], ybar_g[l], mubarsort_g[l], PI.k[l], k0bar_g[l], k1bar_g[l], PI.r[l], h[l], fdiff_g[l])
    critcone_g[l] = outl_g
    hbar_g[l] .= outl_g[:hbar];
  end
  R_d_g = sum(PI.A[l]' * (h[l] .- hbar_g[l]) for l in 1:PI.L);
  R_d_j = sum(PI.A[l]' * (h[l] .- hbar_j[l]) for l in 1:PI.L);
  @assert(isapprox(R_d_g, sum(PI.A[l]' * (h[l] .- hbar_g[l]) for l in 1:1)));
  @assert(isapprox(R_d_j, sum(PI.A[l]' * (h[l] .- hbar_j[l]) for l in 1:1)));
  norm(R_d_g .- R_d_j) #! very different k0bar, k1bar


  # get generalized G
  T = Tf
  m = PI.m[1]
  A = PI.A[1]
  y0 = ytilde[1]
  sig_y = sortperm(y0, rev=true)
  siginv_y = invperm(sig_y)
  y0sort = y0[sig_y]
  ybarsort = similar(y0sort)
  ybar = similar(y0sort)
  r = PI.r[1]
  h = A * d
  k = PI.k[1]

  project_maxksum!(ybarsort, y0sort, r, k, true)
  ybar .= ybarsort[siginv_y]
  musort = get_mu(ybarsort, y0sort, r, k)
  k0bar, k1bar = get_k0k1(ybarsort, k)
  alpha = 1:k0bar
  beta = k0bar+1:k1bar
  gamma = k1bar+1:m
  w = musort[beta]
  beta1 = beta[musort[beta] .== 1]
  beta2 = beta[(musort[beta] .> 0) .* (musort[beta] .< 1)]
  beta3 = beta[musort[beta] .== 0]
  ind_alpha = ones(T, length(alpha))
  A1 = [ind_alpha' w']
  A2 = [zeros(T, length(beta2)-1, length(alpha)) diagm(0 => ones(T, length(beta2)), 1 => -ones(T, length(beta2) - 1))[1:length(beta2)-1,:]]
  B = [A1; A2]
  Gpartial = (I - B' * ((B * B') \ B))
  G = [Gpartial spzeros(T, m - length(gamma), length(gamma)) ; spzeros(T, length(gamma), m - length(gamma)) I]
  Rexplicit = A' * (I - G)[siginv_y,siginv_y] * A


  Cmat(k0, k1) = diagm(0 => ones(k1-k0-1), 1 => -ones(k1-k0-1))[1:end-1,:]
  B_j = Vector{Matrix{Float64}}(undef, PI.L)
  Q_j = Vector{Matrix{Float64}}(undef, PI.L)
  active = zeros(Bool, PI.L)
  for l in 1:PI.L
    if maxksum(ytilde[l], PI.k[l]) > PI.r[l]
      active[l] = 1
      B = [
        ones(k0bar_j[l])'                           mubarsort_j[l][k0bar_j[l]+1:k1bar_j[l]]';
        zeros(k1bar_j[l]-k0bar_j[l]-1, k0bar_j[l])  Cmat(k0bar_j[l], k1bar_j[l])
      ]
      Q = B' * ((B * B') \ B)
      B_j[l] = B
      Q_j[l] = Q
    end
  end
  
  Gpartial = (I - B_j[1]' * ((B_j[1] * B_j[1]') \ B_j[1]))
  alpha = 1:k0bar_j[1]
  beta = k0bar_j[1]+1:k1bar_j[1]
  gamma = k1bar_j[1]+1:PI.m[1]

  nalpha = length(alpha)
  nbeta = length(beta)
  ngamma = length(gamma)
  G = sparse([Gpartial spzeros(PI.m[1] - ngamma, ngamma) ; spzeros(ngamma, PI.m[1] - ngamma) I])
  Rexplicit = PI.A[1]' * (I - G)[siginv_j[1],siginv_j[1]] * PI.A[1]
  # G_j = Vector{SparseMatrixCSC{Float64,Int64}}(undef, PI.L)
  # for l in 1:PI.L
  #   if active[l]
  #     G_j[l] = [
  #       sparse(I(size(Q_j[l],1)) - Q_j[l])                spzeros(size(Q_j[l],1), PI.m[l]-size(Q_j[l],1));
  #       spzeros(size(Q_j[l],1), PI.m[l]-size(Q_j[l],1))'  spdiagm(ones(PI.m[l]-size(Q_j[l],1)))
  #     ]
  #   end
  # end
  # Rexplicit = sigma * sum(PI.A[l]' * (I - G_j[l])[siginv_j[l],siginv_j[l]] * PI.A[l] for l in findall(active))
  # Rexplicit * d

  pp = size(Q_j[1],1)
  Q = sparse( [
    Q_j[1]                   spzeros(pp, PI.m[1]-pp);
    spzeros(PI.m[1]-pp, pp)  spzeros(PI.m[1]-pp, PI.m[1]-pp)
  ])
  siginv = invperm(sig_j[1])
  # T = Q[siginv,siginv] * PI.A[1]
  T = Q * PI.A[1][sig_j[1],:]
  Rexplicit = T'*T
  
  xprev = x;
  Cx = similar(x);
  Ttilde = zeros(sum(PI.m)+PI.L, PI.n);
  Ttilde_rowmax = Base.RefValue{Ti}(0);
  c = [zeros(Rational{Int64}, 3) for l in eachindex(PI.m)];
  a = [zeros(Tf, PI.n) for l in eachindex(PI.m)];
  b = [zeros(Tf, PI.n) for l in eachindex(PI.m)];
  active = Bool[1,0]; # active CVaR constraints
  mks = zeros(Tf, PI.L);
  diagjX = zeros(Tf, PI.n); # diagonal of jacobian of projection onto X
  debug = false
  obj = FunctionValue(Tf, TI.n, TI.objtype); # 2==diagonal
  obj.H .= PI.C;
  phi = FunctionValue(Tf, TI.n, 3); # 2==diagonal
  Hphi!(
    x, xprev, d,
    ytilde, ybar_j, ytmp, sig_j, k0bar_j, k1bar_j,
    utilde, ubar,
    sigma, 0.0,
    phi, obj, PI,
    Cx,
    Ttilde, Ttilde_rowmax, c, a, b, active, diagjX, mks,
    debug
  )
  # R = (phi.H - obj.H .- sigma .* Diagonal(ones(PI.n) .- diagjX)) ./ sigma
  # R_d = R * d
  R = (Ttilde' * Ttilde)
  R_d = (Ttilde' * Ttilde) * d
  norm(R_d - R_d_j)
end


# compare dd and linear dd
hcat(htilde, hbar)
norm(htilde - hbar)


#
# good
#

TI.AP.maxait = 100
TI.SP.maxnit = 500
TI.SP.lin_solver = (:M,1)
n = 5_000
m = [10_000]
k = [100]
CI = ConstraintInstance(12345, n, m, k, TI.scale_row, TI.scale_col, 0.0);
OI = ObjectiveInstance(12345, 2, n);
PI = ProblemInstance(CI, OI);

obj = FunctionValue(Tf, n, 2); # 2==diagonal
obj.H .= PI.C;

x0 = zeros(Tf, n);
lambda0 = [zeros(m[l]) for l in eachindex(m)];
mux0 = zeros(Tf, n);
debug = false;
res = alm(x0, lambda0, mux0, 1.0, 0.0, PI, TI.AP, TI.SP, false)
res_alm = deepcopy(res)

res = solve_cvar_gurobi(PI, options)
calcs_manual!(
  res, obj, PI,
  res[:x], res[:y], res[:z], res[:t],
  res[:lambda], res[:mu]
)
res_gurobi = deepcopy(res)


#
# why maxksum slow
#

inactive = [
  "40_50_596";
  "40_50_667";
  
  "41_20_73";
  "41_20_156";
  
  "41_51_856";
  "41_51_932";
  
  "42_20_499";
  "42_20_579";
  
  "42_46_112";
  "42_46_195";
  
  "43_15_730";
  "43_15_795";
  
  "43_45_535";
  "43_45_614";

  "44_7_832";

  "44_7_904";
  "44_7_983";

  "44_8_160";

  "44_14_559";

  "44_14_751";
]
active = ["44_14_640"]
slow = ["40_52_757", "40_55_708", "40_58_596", "41_1_465", "41_4_372", "41_7_331", "41_10_335", "41_13_389", "41_16_356", "41_18_296", "41_20_304", "41_23_246", "41_26_148", "41_29_0", "41_31_881", "41_34_749", "41_37_772", "41_40_693", "41_43_631", "41_46_568", "41_49_516", "41_50_712", "41_52_78", "41_54_937", "41_57_851", "42_0_711", "42_3_685", "42_6_612", "42_9_569", "42_12_549", "42_15_531", "42_17_880", "42_20_727", "42_23_712", "42_26_627", "42_29_540", "42_32_480", "42_35_446", "42_38_463", "42_41_140", "42_46_342", "42_49_274", "42_52_200", "42_55_73", "42_57_984", "43_0_869", "43_3_801", "43_6_853", "43_9_824", "43_15_938", "43_18_799", "43_21_718", "43_24_583", "43_27_594", "43_30_508", "43_33_449", "43_36_416", "43_39_370", "43_41_653", "43_45_762", "43_48_662", "43_51_700", "43_54_655", "43_57_656", "44_0_193"] # among others
med = ["42_43_658", "43_14_884", "43_15_317", "43_43_896", "43_44_790", "44_2_699", "44_4_308", "44_5_735", "44_6_468", "44_8_315", "44_10_376", "44_12_310", "44_13_99", "44_14_903", "44_15_473"]
fast = ["42_20_225", "42_20_419", "44_7_82", "44_7_402", "44_7_599", "44_7_832", "44_8_63", "44_13_839", "44_14_136", "44_14_327", "44_14_640"]

id = slow[1]
ytildesort_slow = vec(readdlm("slowmaxksum/ytildesort_$(id).csv"));
ytildesortbar_slow = similar(ytildesort_slow)
r = vec(readdlm("slowmaxksum/r_$(id).csv"))[1];
k = Int(vec(readdlm("slowmaxksum/k_$(id).csv"))[1]);
@views out = project_maxksum!(ytildesortbar_slow, ytildesort_slow, r, k);
get_k0k1!(ytildesortbar_slow, k)
(k0s, k1s, thetas, lambdas, rhos) = out[2]
n = 131072
writedlm("k0s.csv", k0s[1:n-k])
writedlm("k1s.csv", k1s[1:n-k])
writedlm("thetas.csv", thetas[1:n-k])
writedlm("lambdas.csv", lambdas[1:n-k])
writedlm("rhos.csv", rhos[1:n-k])

k0s = Int.(vec(readdlm(PROJPATH * "k0s.csv")))
k1s = Int.(vec(readdlm(PROJPATH * "k1s.csv")))
thetas = vec(readdlm(PROJPATH * "thetas.csv"))
lambdas = vec(readdlm(PROJPATH * "lambdas.csv"))
rhos = vec(readdlm(PROJPATH * "rhos.csv"))
ax.cla(); ax.plot(k1s, lambdas, label="lam"); ax.legend()
ax.cla(); ax.plot(k1s, thetas, label="theta"); ax.legend()


#
# slow ssn convergence
#

rt = readdlm("trace.txt", '│')
gphi = [parse(Float64,x[1]) for x in split.(rt[3:end-1,10],"│")]
ax.cla(); ax.semilogy(collect(1:length(gphi)),gphi, label="||gphi||"); ax.legend(); ax.set_xlabel("iter")
=#