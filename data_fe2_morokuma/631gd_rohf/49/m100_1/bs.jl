using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf

@load "data_cmf_49_new.jld2"

M = 50


init_fspace =  [(5, 0), (4, 4), (0, 5), (12, 12), (0, 0)]
init_fspace = FermiCG.FockConfig(init_fspace)
# cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,1,1], init_fspace, max_roots=M, verbose=1);
@load "cluster_bases_49.jld2" 
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=6
# start by defining P/Q spaces
p_spaces = Vector{ClusterSubspace}()

    
ssi = ClusterSubspace(clusters[1])
add_subspace!(ssi, (5,0), 1:1)
add_subspace!(ssi, (4,1), 1:1)
add_subspace!(ssi, (3,2), 1:1)
add_subspace!(ssi, (2,3), 1:1)
add_subspace!(ssi, (1,4), 1:1)
add_subspace!(ssi, (0,5), 1:1)
push!(p_spaces, ssi)

ssi = ClusterSubspace(clusters[2])
add_subspace!(ssi, init_fspace[2], 1:1)
push!(p_spaces, ssi)

ssi = ClusterSubspace(clusters[3])
add_subspace!(ssi, (5,0), 1:1)
add_subspace!(ssi, (4,1), 1:1)
add_subspace!(ssi, (3,2), 1:1)
add_subspace!(ssi, (2,3), 1:1)
add_subspace!(ssi, (1,4), 1:1)
add_subspace!(ssi, (0,5), 1:1)
push!(p_spaces, ssi)

ssi = ClusterSubspace(clusters[4])
add_subspace!(ssi, init_fspace[4], 1:1)
push!(p_spaces, ssi)

ssi = ClusterSubspace(clusters[5])
add_subspace!(ssi, init_fspace[5], 1:1)
push!(p_spaces, ssi)


v = FermiCG.BSstate(clusters,p_spaces, cluster_bases, R=6) 
na = sum([i[1] for i in init_fspace]) 
nb = sum([i[2] for i in init_fspace]) 

FermiCG.fill_p_space!(v, na, nb)
# v = FermiCG.BSstate(clusters, init_fspace, cluster_bases, R=6)
FermiCG.eye!(v)
display(v)
e_ci, v_ci = FermiCG.ci_solve(v, cluster_ops, clustered_ham, solver="davidson",max_iter = 100);

v_bst = FermiCG.BSTstate(v_ci, thresh=1e-3)

display(v_bst)
FermiCG.randomize!(v_bst)
FermiCG.orthonormalize!(v_bst)

display(v_bst)
σ = FermiCG.build_compressed_1st_order_state(v_bst, cluster_ops, clustered_ham, 
                                    nbody=4,
                                    thresh=1e-3)
σ = FermiCG.compress(σ, thresh=1e-3)
v2 = BSTstate(σ,R=6)
FermiCG.eye!(v2)
e_ci, v2 = FermiCG.ci_solve(v2, cluster_ops, clustered_ham);
e_var, v_var = FermiCG.block_sparse_tucker(v2, cluster_ops, clustered_ham,
                                               max_iter    = 20,
                                               nbody       = 4,
                                               H0          = "Hcmf",
                                               thresh_var  = 1e-2,
                                               thresh_foi  = 1e-4,
                                               thresh_pt   = 1e-3,
                                               ci_conv     = 1e-5,
                                               do_pt       = true,
                                               resolve_ss  = false,
                                               tol_tucker  = 1e-4,
                                               solver      = "davidson")
