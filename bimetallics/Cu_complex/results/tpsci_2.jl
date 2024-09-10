using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using LinearAlgebra

@load "data_cmf_32_Cu_dimer.jld2"

M = 100
init_fspace =  [(5, 4),  (3, 3), (4, 5), (3, 3)]


init_fspace =  FockConfig(init_fspace)
@load "tpsci_out.jld2"
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3], init_fspace, max_roots=M, verbose=1);

#clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=2
rotations = FermiCG.hosvd(v0a, cluster_ops)
for ci in clusters
     FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
     FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 1e-5,
                                thresh_spin  = 1e-5,
                                thresh_foi   = 1e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
@save "/projects/nmayhall_lab/arnab/dmrg/Cu_complex/tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 
rotations = FermiCG.hosvd(v0a, cluster_ops)
for ci in clusters
     FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
     FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 9e-6,
                                thresh_spin  = 9e-6,
                                thresh_foi   = 1e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
@save "/projects/nmayhall_lab/arnab/dmrg/Cu_complex/tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 
rotations = FermiCG.hosvd(v0a, cluster_ops)
for ci in clusters
     FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
     FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 8e-6,
                                thresh_spin  = 8e-6,
                                thresh_foi   = 1e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
@save "/projects/nmayhall_lab/arnab/dmrg/Cu_complex/tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 
rotations = FermiCG.hosvd(v0a, cluster_ops)
for ci in clusters
     FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
     FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 7e-6,
                                thresh_spin  = 7e-6,
                                thresh_foi   = 1e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
@save "/projects/nmayhall_lab/arnab/dmrg/Cu_complex/tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 
rotations = FermiCG.hosvd(v0a, cluster_ops)
for ci in clusters
     FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
     FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 6e-6,
                                thresh_spin  = 6e-6,
                                thresh_foi   = 1e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
@save "/projects/nmayhall_lab/arnab/dmrg/Cu_complex/tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 
rotations = FermiCG.hosvd(v0a, cluster_ops)
for ci in clusters
     FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
     FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 5e-6,
                                thresh_spin  = 5e-6,
                                thresh_foi   = 1e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
@save "/projects/nmayhall_lab/arnab/dmrg/Cu_complex/tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 
