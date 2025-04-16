using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using LinearAlgebra

@load "data_cmf_M10_Fe_dimer.jld2"

M = 30
init_fspace =  [(5, 0), (0,5), (3, 3), (3, 3), (3, 3)]


init_fspace =  FockConfig(init_fspace)

cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3,3], init_fspace, max_roots=M, verbose=1);

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);



nroots=10
ci_vector1 = FermiCG.TPSCIstate(clusters, FockConfig([(5, 1), (0,5), (3, 3), (3, 3), (3, 3)]), R=nroots)
ci_vector2 = FermiCG.TPSCIstate(clusters, FockConfig([(0, 5), (5,1), (3, 3), (3, 3), (3, 3)]), R=nroots)

ci_vector1 = FermiCG.add_spin_focksectors(ci_vector1)
ci_vector2 = FermiCG.add_spin_focksectors(ci_vector2)

ci_vector = ci_vector1 + ci_vector2

# ci_vector3 = FermiCG.TPSCIstate(clusters, FockConfig([(4,2), (0,5), (3, 3), (3, 3), (3, 3)]), R=nroots)
# ci_vector4 = FermiCG.TPSCIstate(clusters, FockConfig([(3,3), (0,5), (3, 3), (3, 3), (3, 3)]), R=nroots)
# ci_vector5 = FermiCG.TPSCIstate(clusters, FockConfig([(5,1), (1,4), (3, 3), (3, 3), (3, 3)]), R=nroots)
# ci_vector6 = FermiCG.TPSCIstate(clusters, FockConfig([(5,1), (2,3), (3, 3), (3, 3), (3, 3)]), R=nroots)

# ci_vector7 = FermiCG.TPSCIstate(clusters, FockConfig([(0,5),(4,2),  (3, 3), (3, 3), (3, 3)]), R=nroots)
# ci_vector8 = FermiCG.TPSCIstate(clusters, FockConfig([(0,5),(3,3),  (3, 3), (3, 3), (3, 3)]), R=nroots)
# ci_vector9 = FermiCG.TPSCIstate(clusters, FockConfig([(1,4),(5,1),  (3, 3), (3, 3), (3, 3)]), R=nroots)
# ci_vector10 = FermiCG.TPSCIstate(clusters, FockConfig([(2,3),(5,1),  (3, 3), (3, 3), (3, 3)]), R=nroots)

# ci_vector__= ci_vector3 +ci_vector4+ ci_vector5  + ci_vector9 + ci_vector10
# ci_vector__ = FermiCG.add_spin_focksectors(ci_vector__)
# ci_vector= ci_vector_ + ci_vector__
display(ci_vector)

eci, v = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);

e0a, v0a = FermiCG.tpsci_ci(v, cluster_ops, clustered_ham,
                                thresh_cipsi = 6e-3,
                                thresh_foi   = 5e-5,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);
@save "tpsci_out.jld2"  cluster_bases clustered_ham init_fspace  cluster_ops

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 1e-3,
                                thresh_spin  = 1e-3,
                                thresh_foi   = 5e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

@save "tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 1e-3,
                                thresh_spin  = 1e-3,
                                thresh_foi   = 1e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

@save "tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 9e-4,
                                thresh_spin  = 9e-4,
                                thresh_foi   = 1e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

@save "tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 8e-4,
                                thresh_spin  = 8e-4,
                                thresh_foi   = 1e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

@save "tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 7e-4,
                                thresh_spin  = 7e-4,
                                thresh_foi   = 1e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

@save "tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
