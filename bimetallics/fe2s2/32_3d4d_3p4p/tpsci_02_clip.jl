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

M = 200 

nroots=6
init_fspace =  init_fspace =  [(5, 0), (3, 3), (3, 3), (0, 5)]
init_fspace = FermiCG.FockConfig(init_fspace)
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3], init_fspace, max_roots=M, verbose=1);
@load "data_tpsci_200.jld2"
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);
#FermiCG.clip!(v0a,thresh=5e-5)
#e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
#print_pt2(e2a, nroots)
FermiCG.clip!(v0a,thresh=6e-5)
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)
FermiCG.clip!(v0a,thresh=7e-5)
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)
FermiCG.clip!(v0a,thresh=8e-5)
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)
