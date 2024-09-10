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

M = 150 


init_fspace =  [(5, 0), (4, 4), (0, 5), (12, 12), (0, 0)]
init_fspace = FermiCG.FockConfig(init_fspace)
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,1,1], init_fspace, max_roots=M, verbose=1);
@load "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/cluster_bases_150.jld2" 
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=6


ci_vector = FermiCG.TPSCIstate(clusters, init_fspace, R=nroots)

ci_vector = FermiCG.add_spin_focksectors(ci_vector)

eci, v = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);
@load "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/data_tpsci_150.jld2"
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 3.5e-4, 
                            thresh_spin  = 5e-4,
                            thresh_foi   = 1e-5,
                            max_mem_ci   = 900);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/data_tpsci_150.jld2" clusters init_fspace ints cluster_bases v0a
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 3.2e-4, 
                            thresh_spin  = 4.5e-4,
                            thresh_foi   = 1e-5,
                            max_mem_ci   = 900);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)
@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/data_tpsci_150.jld2" clusters init_fspace ints cluster_bases v0a
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 3e-4, 
                            thresh_spin  = 4e-4,
                            thresh_foi   = 1e-5,
                            max_mem_ci   = 900);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)


@save "/home/arnab22/FermiCG-data/bimetallics/fe2_morokuma/rohf.631gd/49/data_tpsci_150.jld2" clusters init_fspace ints cluster_bases v0a


