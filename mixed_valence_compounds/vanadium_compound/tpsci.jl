using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using LinearAlgebra

@load "data_cmf_M6_V_dimer.jld2"

M = 30
init_fspace =  [(3, 0), (0,3)]


init_fspace =  FockConfig(init_fspace)

cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [5,5], init_fspace, max_roots=M, verbose=1);
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=18

ci_vector1= FermiCG.TPSCIstate(clusters, FockConfig([(2, 0), (0,3)]), R=nroots)
ci_vector2= FermiCG.TPSCIstate(clusters, FockConfig([(0, 3), (2,0)]), R=nroots)
ci_vector1 = FermiCG.add_spin_focksectors(ci_vector1)
ci_vector2 = FermiCG.add_spin_focksectors(ci_vector2)


# ci_vector3= FermiCG.TPSCIstate(clusters, FockConfig([(1, 1), (0,3)]), R=nroots)
# ci_vector4= FermiCG.TPSCIstate(clusters, FockConfig([(1, 1), (1,2)]), R=nroots)
# ci_vector5= FermiCG.TPSCIstate(clusters, FockConfig([(1, 1), (2,1)]), R=nroots)
# ci_vector6= FermiCG.TPSCIstate(clusters, FockConfig([(1, 1), (3,0)]), R=nroots)


# ci_vector7= FermiCG.TPSCIstate(clusters, FockConfig([(2, 0), (1,2)]), R=nroots)
# ci_vector8= FermiCG.TPSCIstate(clusters, FockConfig([(2, 0), (2,1)]), R=nroots)
# ci_vector9= FermiCG.TPSCIstate(clusters, FockConfig([(0, 2), (1,2)]), R=nroots)
# ci_vector10= FermiCG.TPSCIstate(clusters, FockConfig([(0, 2), (2,1)]), R=nroots)

# ci_vector11= FermiCG.TPSCIstate(clusters, FockConfig([(0, 3), (1,1)]), R=nroots)
# ci_vector12= FermiCG.TPSCIstate(clusters, FockConfig([(1, 2), (1,1)]), R=nroots)
# ci_vector13= FermiCG.TPSCIstate(clusters, FockConfig([(2, 1), (1,1)]), R=nroots)
# ci_vector14= FermiCG.TPSCIstate(clusters, FockConfig([(3, 0), (1,1)]), R=nroots)
# ci_vector15= FermiCG.TPSCIstate(clusters, FockConfig([(1, 2), (2,0)]), R=nroots)
# ci_vector16= FermiCG.TPSCIstate(clusters, FockConfig([(2, 1), (2,0)]), R=nroots)
# ci_vector17= FermiCG.TPSCIstate(clusters, FockConfig([(1, 2), (0,2)]), R=nroots)
# ci_vector18= FermiCG.TPSCIstate(clusters, FockConfig([(2, 1), (0,2)]), R=nroots)

ci_vector = ci_vector1 + ci_vector2#+ ci_vector3 + ci_vector4+ ci_vector5 + ci_vector6+ ci_vector7 + ci_vector8+ ci_vector9 + ci_vector10+ ci_vector11 + ci_vector12+ ci_vector13 + ci_vector14+ ci_vector15 + ci_vector16+ ci_vector17 + ci_vector18
display(ci_vector) 


eci, v = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);
e0a, v0a = FermiCG.tpsci_ci(v, cluster_ops, clustered_ham,
                                thresh_cipsi = 5e-3,
                                thresh_spin  = 4e-3,
                                thresh_foi   = 1e-5,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

@save "tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 

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
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham,
                                thresh_cipsi = 1e-5,
                                thresh_spin  = 1e-5,
                                thresh_foi   = 1e-6,
                                ci_max_iter  = 150,
                                conv_thresh  = 1e-6,
                                nbody        = 4,
                                max_mem_ci   = 200);

@save "tpsci_out.jld2" ci_vector v0a e0a cluster_bases clustered_ham init_fspace 

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-8)
