using QCBase
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf

function print_pt2(e, nroots)
    @printf("*E(PT2): ")
    for i in e
        @printf("%14.12f ", i)
    end
    println()
end

@load "data_cmf.jld2"
@load "data_tpsci.jld2"
M = 100 


init_fspace =  [(5, 0), (4, 4), (0, 5), (12, 12), (0, 0)]
init_fspace = FermiCG.FockConfig(init_fspace)
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,1,1], init_fspace, max_roots=M, verbose=1);

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=6



e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 3e-4, 
                            thresh_spin  = 2e-4,
                            thresh_foi   = 1e-5,
                            max_mem_ci   = 800);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1.0e-8)
print_pt2(e2a, nroots)
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/data_tpsci.jld2" clusters init_fspace ints cluster_bases v0a
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 2.5e-4, 
                            thresh_spin  = 1e-4,
                            thresh_foi   = 1e-5,
                            max_mem_ci   = 800);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham,thresh_foi=1.0e-8)
print_pt2(e2a, nroots)
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/data_tpsci.jld2" clusters init_fspace ints cluster_bases v0a
