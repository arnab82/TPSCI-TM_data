using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf
using ActiveSpaceSolvers
C = npzread("/home/arnabbachhar/workspace/FermiCG-data/bimetallics/fe2_morokuma/26_3d4d_2p3p_3d4d/mo_coeffs_mb.npy")
h0 = npzread("/home/arnabbachhar/workspace/FermiCG-data/bimetallics/fe2_morokuma/26_3d4d_2p3p_3d4d/ints_h0_mb.npy")
h1 = npzread("/home/arnabbachhar/workspace/FermiCG-data/bimetallics/fe2_morokuma/26_3d4d_2p3p_3d4d/ints_h1_mb.npy")
h2 = npzread("/home/arnabbachhar/workspace/FermiCG-data/bimetallics/fe2_morokuma/26_3d4d_2p3p_3d4d/ints_h2_mb.npy")
ints = InCoreInts(h0, h1, h2)

Pa = npzread("/home/arnabbachhar/workspace/FermiCG-data/bimetallics/fe2_morokuma/26_3d4d_2p3p_3d4d/Pa_mb.npy")
Pb = npzread("/home/arnabbachhar/workspace/FermiCG-data/bimetallics/fe2_morokuma/26_3d4d_2p3p_3d4d/Pb_mb.npy")
@printf(" Input energy:    %12.8f\n", compute_energy(ints, RDM1(Pa, Pb)))


init_fspace=  [(5, 1), (3, 1), (5, 1)]
clusters   =  [[1, 2, 3, 4, 5], [6, 7, 8], [9, 10, 11, 12, 13]]



clusters = [MOCluster(i, collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

d1 = RDM1(n_orb(ints))

# ansatze = [FCIAnsatz(5,5,0), FCIAnsatz(3,3,3),FCIAnsatz(5,5,0)]
e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace,d1,
                           maxiter_oo   = 700, 
                           maxiter_ci   = 200, 
                           maxiter_d1   = 200, 
                           verbose      = 0, 
                           tol_oo       = 1e-6, 
                           tol_d1       = 1e-9, 
                           tol_ci       = 1e-11, 
                          sequential   = true, 
                           alpha        = .1,
                           diis_start   = 1,
                           max_ss_size  = 36)

ints = orbital_rotation(ints, U)
C = C*U

npzwrite("Ccmf_fe2_13_ano-rcc-mb.npy", C)

@save "data_cmf_fe2_ano-rcc-mb.jld2" clusters init_fspace ints d1 e_cmf U 