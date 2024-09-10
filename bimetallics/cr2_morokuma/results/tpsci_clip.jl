using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


@load "data_cmf_29_cr2.jld2"

M = 250

init_fspace = FockConfig([(6,3), (4, 4), (3, 6)])

#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3], init_fspace, max_roots=M, verbose=1);
@load "data_tpsci_6e5.jld2"
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);
@time e2 = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=4);
nroots=4
FermiCG.clip!(v0a, thresh=0.00007)
e0, v0 = FermiCG.tps_ci_direct(v0a, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v0, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=4);

clustered_S2 = FermiCG.extract_S2(v0.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v0, cluster_ops, clustered_S2)
FermiCG.clip!(v0a, thresh=0.00008)
e0, v01 = FermiCG.tps_ci_direct(v0a, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v01, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=4);

clustered_S2 = FermiCG.extract_S2(v01.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v01, cluster_ops, clustered_S2)
FermiCG.clip!(v0a, thresh=0.00009)
e0, v02 = FermiCG.tps_ci_direct(v0a, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v02, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=4);

clustered_S2 = FermiCG.extract_S2(v02.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v02, cluster_ops, clustered_S2)
FermiCG.clip!(v0a, thresh=0.0001)
e0, v03 = FermiCG.tps_ci_direct(v0a, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v03, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=4);

clustered_S2 = FermiCG.extract_S2(v03.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v03, cluster_ops, clustered_S2)

FermiCG.clip!(v0a, thresh=0.0002)
e0, v04 = FermiCG.tps_ci_direct(v0a, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v04, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=4);

clustered_S2 = FermiCG.extract_S2(v04.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v04, cluster_ops, clustered_S2)
FermiCG.clip!(v0a, thresh=0.0003)
e0, v11 = FermiCG.tps_ci_direct(v0a, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v11, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=4);

clustered_S2 = FermiCG.extract_S2(v11.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v11, cluster_ops, clustered_S2)

FermiCG.clip!(v0a, thresh=0.0004)
e0, v11_ = FermiCG.tps_ci_direct(v0a, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v11_, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=4);

clustered_S2 = FermiCG.extract_S2(v11_.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v11_, cluster_ops, clustered_S2)
FermiCG.clip!(v0a, thresh=0.0005)
e0, v12 = FermiCG.tps_ci_direct(v0a, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v12, cluster_ops, clustered_ham, thresh_foi=1e-7,verbose=4);
