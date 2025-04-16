using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using LinearAlgebra

@load "data_cmf_32_Mn_Cy.jld2"

M = 30
init_fspace =  [(5, 4), (0, 5), (3, 3)]
init_fspace =  FockConfig(init_fspace)

cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3], init_fspace, max_roots=M, verbose=1);
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops   = FermiCG.compute_cluster_ops(cluster_bases, ints);
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=6
ci_vector = FermiCG.TPSCIstate(clusters, init_fspace, R=nroots)

FermiCG.add_fockconfig!(ci_vector, FockConfig([(5, 4), (1, 4), (3, 3)]))
FermiCG.add_fockconfig!(ci_vector, FockConfig([(5, 4), (2, 3), (3, 3)]))
ci_vector = FermiCG.add_spin_focksectors(ci_vector)
display(ci_vector)
eci, v    = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);

# @load "tpsci_out.jld2"

e0a, v0a  = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 1e-3,
                                thresh_spin  = 1e-3,
                                thresh_foi   = 5e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

@save "tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 

# e2a       = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
e0a, v0a  = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 6e-4,
                                thresh_spin  = 4e-4,
                                thresh_foi   = 5e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

@save "tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 

e2a       = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
e0a, v0a  = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 5e-4,
                                thresh_spin  =3e-4,
                                thresh_foi   = 4e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

@save "tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 

e2a       = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
e0a, v0a  = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 4e-4,
                                thresh_spin  = 2e-4,
                                thresh_foi   = 3e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

@save "tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 

e2a       = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
e0a, v0a  = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 3e-4,
                                thresh_spin  = 1e-4,
                                thresh_foi   = 2e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

@save "tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 

e2a       = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)