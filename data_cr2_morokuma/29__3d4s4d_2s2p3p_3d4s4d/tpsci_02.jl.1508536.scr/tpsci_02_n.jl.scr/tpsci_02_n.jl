using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


@load "/home/arnab22/FermiCG-data/bimetallics/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/data_cmf_29_cr2.jld2"

M = 200 

init_fspace = FockConfig([(6,3), (4, 4), (3, 6)])

#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3], init_fspace, max_roots=M, verbose=1);
@load "/home/arnab22/FermiCG-data/bimetallics/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/tpsci_02.jl.1508536.scr/cluster_bases_02.jld2" 
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=4


ci_vector = FermiCG.TPSCIstate(clusters, init_fspace, R=nroots)

ci_vector = FermiCG.add_spin_focksectors(ci_vector)

eci, v = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);

e0a, v0a = FermiCG.tpsci_ci(v, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 4e-4, 
                            thresh_foi   = 1e-6,
                            thresh_asci  = -1,
			    max_mem_ci   = 200.0);

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
@save "/home/arnab22/FermiCG-data/bimetallics/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/tpsci_02.jl.1508536.scr/data_tpsci_02.jld2" e0a v0a e2a
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 2e-4,
                            thresh_foi   = 1e-6,
                            thresh_asci  = -1,
			    max_mem_ci   = 200.0);

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
@save "/home/arnab22/FermiCG-data/bimetallics/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/tpsci_02.jl.1508536.scr/data_tpsci_02.jld2" e0a v0a e2a
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 1e-4,
                            thresh_foi   = 1e-6,
                            thresh_asci  = -1,
			    max_mem_ci   = 200.0);

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
@save "/home/arnab22/FermiCG-data/bimetallics/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/tpsci_02.jl.1508536.scr/data_tpsci_02.jld2" e0a v0a e2a
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 8e-5,
                            thresh_foi   = 1e-6,
                            thresh_asci  = -1,
                            max_mem_ci   = 200.0);

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
@save "/home/arnab22/FermiCG-data/bimetallics/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/tpsci_02.jl.1508536.scr/data_tpsci_02.jld2" e0a v0a e2a
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 6e-5,
                            thresh_foi   = 1e-6,
                            thresh_asci  = -1,
                            max_mem_ci   = 200.0);

e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
@save "/home/arnab22/FermiCG-data/bimetallics/cr2_morokuma/29__3d4s4d_2s2p3p_3d4s4d/tpsci_02.jl.1508536.scr/data_tpsci_02.jld2" e0a v0a e2a
