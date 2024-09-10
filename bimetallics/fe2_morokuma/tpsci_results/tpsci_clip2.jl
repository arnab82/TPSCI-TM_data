using FermiCG
using JLD2
using LinearAlgebra
using Printf
using QCBase
using RDM

@load  "data_cmf.jld2"
ref_fspace = FockConfig(init_fspace)
ecore = ints.h0
@load  "data_tpsci.jld2"
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);



e0, v11 = FermiCG.tps_ci_direct(v0a, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v11, cluster_ops, clustered_ham, thresh_foi=1e-8);
clustered_S2 = FermiCG.extract_S2(v11.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v11, cluster_ops, clustered_S2)
@printf(" %5s %12.8f %12.8f %12.8f\n",1, e0[1]+ecore,e2[1]+ecore, abs(s2[1]))
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/clip/tucker_clip_0.0004.jld2" clusters cluster_bases e0 v11 e2 ecore s2

FermiCG.clip!(v11, thresh=0.0006)
e0, v12 = FermiCG.tps_ci_direct(v11, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v12, cluster_ops, clustered_ham, thresh_foi=1e-8);
clustered_S2 = FermiCG.extract_S2(v12.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v12, cluster_ops, clustered_S2)
@printf(" %5s %12.8f %12.8f %12.8f\n",1, e0[1]+ecore,e2[1]+ecore, abs(s2[1]))
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/clip/tucker_clip_0.0006.jld2" clusters cluster_bases e0 v12 e2 ecore s2

FermiCG.clip!(v12, thresh=0.0008)
e0, v13 = FermiCG.tps_ci_direct(v12, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v13, cluster_ops, clustered_ham, thresh_foi=1e-8);
clustered_S2 = FermiCG.extract_S2(v13.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v13, cluster_ops, clustered_S2)
@printf(" %5s %12.8f %12.8f %12.8f\n",1, e0[1]+ecore,e2[1]+ecore, abs(s2[1]))
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/clip/tucker_clip_0.0008.jld2" clusters cluster_bases e0 v13 e2 ecore s2

FermiCG.clip!(v13, thresh=0.001)
e0, v14 = FermiCG.tps_ci_direct(v13, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v14, cluster_ops, clustered_ham, thresh_foi=1e-8);
clustered_S2 = FermiCG.extract_S2(v14.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v14, cluster_ops, clustered_S2)
@printf(" %5s %12.8f %12.8f %12.8f\n",1, e0[1]+ecore,e2[1]+ecore, abs(s2[1]))
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/clip/tucker_clip_0.001.jld2" clusters cluster_bases e0 v14 e2 ecore s2

FermiCG.clip!(v14, thresh=0.002)
e0, v15 = FermiCG.tps_ci_direct(v14, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v15, cluster_ops, clustered_ham, thresh_foi=1e-8);
clustered_S2 = FermiCG.extract_S2(v15.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v15, cluster_ops, clustered_S2)
@printf(" %5s %12.8f %12.8f %12.8f\n",1, e0[1]+ecore,e2[1]+ecore, abs(s2[1]))
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/clip/tucker_clip_0.002.jld2" clusters cluster_bases e0 v15 e2 ecore s2

FermiCG.clip!(v15, thresh=0.004)
e0, v16 = FermiCG.tps_ci_direct(v15, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v16, cluster_ops, clustered_ham, thresh_foi=1e-8);
clustered_S2 = FermiCG.extract_S2(v16.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v16, cluster_ops, clustered_S2)
@printf(" %5s %12.8f %12.8f %12.8f\n",1, e0[1]+ecore,e2[1]+ecore, abs(s2[1]))
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/clip/tucker_clip_0.004.jld2" clusters cluster_bases e0 v16 e2 ecore s2

