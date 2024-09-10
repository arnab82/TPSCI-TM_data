using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


@load "/home/arnab22/FermiCG-data/bimetallics/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/data_cmf_29_cr2.jld2"

M = 250 

init_fspace = FockConfig([(6,3), (4, 4), (3, 6)])

#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3], init_fspace, max_roots=M, verbose=1);
@load "data_bst_7e3_250.jld2"
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=4


ci_vector = FermiCG.TPSCIstate(clusters, init_fspace, R=nroots)

ci_vector = FermiCG.add_spin_focksectors(ci_vector)

eci, v = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);


e0a, v0a = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 3e-4,
                            thresh_foi   = 1e-6,
                            thresh_asci  = -1,
                            max_mem_ci   = 200.0);

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-7,verbose=4)
@save "/home/arnab22/FermiCG-data/bimetallics/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/data_tpsci.jld2" cluster_bases e0a v0a e2a
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 2e-4, 
                            thresh_foi   = 1e-6,
                            thresh_asci  = -1,
			    max_mem_ci   = 200.0);

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1e-7,verbose=4)
@save "/home/arnab22/FermiCG-data/bimetallics/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/data_tpsci.jld2" cluster_bases e0a v0a e2a 
v11=FermiCG.clip!(v0a, thresh=0.0003)
e0, v11 = FermiCG.tps_ci_direct(v11, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v11, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=4);

clustered_S2 = FermiCG.extract_S2(v11.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v11, cluster_ops, clustered_S2)

v11_=FermiCG.clip!(v0a, thresh=0.0004)
e0, v11_ = FermiCG.tps_ci_direct(v11_, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v11_, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=4);

clustered_S2 = FermiCG.extract_S2(v11_.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v11_, cluster_ops, clustered_S2)
v12=FermiCG.clip!(v0a, thresh=0.0005)
e0, v12 = FermiCG.tps_ci_direct(v12, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v12, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=4);

clustered_S2 = FermiCG.extract_S2(v12.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v12, cluster_ops, clustered_S2)
v13=FermiCG.clip!(v0a, thresh=0.0006)
e0, v13 = FermiCG.tps_ci_direct(v13, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v13, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=4);

clustered_S2 = FermiCG.extract_S2(v13.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v13, cluster_ops, clustered_S2)
