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

M = 100 


init_fspace =  init_fspace =  [(5, 0), (3, 3), (3, 3), (0, 5)]
init_fspace = FermiCG.FockConfig(init_fspace)
#cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3], init_fspace, max_roots=M, verbose=1);
@load "data_tpsci.jld2"
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=6


ci_vector = FermiCG.TPSCIstate(clusters, init_fspace, R=nroots)

ci_vector = FermiCG.add_spin_focksectors(ci_vector)

eci, v = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);



e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 1e-4, 
                            thresh_spin  = 1e-4,
                            thresh_foi   = 1e-6,
                            max_iter        = 20,
                            conv_thresh     = 1e-6,
                            max_mem_ci   = 200);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)
@save "/projects/nmayhall_lab/arnab/dmrg/fe2s2/32__3d4d_3p4p/data_tpsci.jld2" clusters init_fspace  cluster_bases v0a e2a e0a
rotations = FermiCG.hosvd(v0a, cluster_ops)
for ci in clusters
     FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
     FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 9e-5,
                            thresh_spin  = 9e-5,
                            thresh_foi   = 1e-6,
                            max_iter        = 20,
                            conv_thresh     = 1e-6,
                            max_mem_ci   = 300);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)
@save "/projects/nmayhall_lab/arnab/dmrg/fe2s2/32__3d4d_3p4p/data_tpsci.jld2" clusters init_fspace  cluster_bases v0a e2a e0a
rotations = FermiCG.hosvd(v0a, cluster_ops)
for ci in clusters
     FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
     FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 8e-5,
                            thresh_spin  = 8e-5,
                            thresh_foi   = 1e-6,
                            max_iter        = 20,
                            conv_thresh     = 1e-6,
                            max_mem_ci   = 200);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)
@save "/projects/nmayhall_lab/arnab/dmrg/fe2s2/32__3d4d_3p4p/data_tpsci.jld2" clusters init_fspace  cluster_bases v0a e2a e0a
rotations = FermiCG.hosvd(v0a, cluster_ops)
for ci in clusters
     FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
     FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 7e-5,
                            thresh_spin  = 7e-5,
                            thresh_foi   = 1e-6,
                            max_iter        = 20,
                            conv_thresh     = 1e-6,
                            max_mem_ci   = 200);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)
@save "/projects/nmayhall_lab/arnab/dmrg/fe2s2/32__3d4d_3p4p/data_tpsci.jld2" clusters init_fspace  cluster_bases v0a e2a e0a
rotations = FermiCG.hosvd(v0a, cluster_ops)
for ci in clusters
     FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
     FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end
e0a, v0a = FermiCG.tpsci_ci(v0a, cluster_ops, clustered_ham, incremental=true,
                            thresh_cipsi = 6e-5,
                            thresh_spin  = 6e-5,
                            thresh_foi   = 1e-6,
                            max_iter        = 20,
                            conv_thresh     = 1e-6,
                            max_mem_ci   = 200);
e2a = FermiCG.compute_pt2_energy(v0a, cluster_ops, clustered_ham)
print_pt2(e2a, nroots)
@save "/projects/nmayhall_lab/arnab/dmrg/fe2s2/32__3d4d_3p4p/data_tpsci.jld2" clusters init_fspace  cluster_bases v0a e2a e0a
