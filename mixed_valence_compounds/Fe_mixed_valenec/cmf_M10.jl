using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf
using ActiveSpaceSolvers

C = npzread("mo_coeffs_M10.npy")
h0 = npzread("ints_h0_M10.npy")
h1 = npzread("ints_h1_M10.npy")
h2 = npzread("ints_h2_M10.npy")
ints = InCoreInts(h0, h1, h2)

Pa = npzread("Pa_M10.npy")
Pb = npzread("Pb_M10.npy")
@printf(" Input energy:    %12.8f\n", compute_energy(ints, RDM1(Pa, Pb)))


init_fspace =  [(5, 0), (5, 0), (3, 3), (3, 3), (3, 3)]
clusters    =  [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22], [23, 24, 25, 26, 27, 28, 29], [30, 31, 32, 33, 34, 35, 36], [37, 38, 39, 40, 41, 42, 43]]
clusters = [MOCluster(i, collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

d1 = RDM1(n_orb(ints))


# # Do CMF
ansatze = [FCIAnsatz(11, 5,0), FCIAnsatz(11,5,0),FCIAnsatz(7,3,3),FCIAnsatz(7,3,3),FCIAnsatz(7,3,3)]

# @time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace,ansatze,d1,max_iter_oo=20,
#                            gconv=1e-6, 
#                            tol_d1=1e-9, 
#                            tol_ci=1e-11,
#                            verbose=0, 
#                            sequential=true)
# ints = orbital_rotation(ints, U)
# C= C*U
# d1=orbital_rotation(d1, U)
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton(ints, clusters, init_fspace,ansatze,d1, maxiter_oo = 400,
                           tol_oo=1e-6, 
                           tol_d1=1e-9, 
                           tol_ci=1e-11,
                           step_trust_region=2.5,
                           verbose=0, 
                           trust_region=true,
                           sequential=true)
ints = orbital_rotation(ints, U)
C= C*U
d1=orbital_rotation(d1, U)
npzwrite("Ccmf_M10.npy", C)
@save "data_cmf_M10_Fe_dimer.jld2" clusters init_fspace ints d1 e_cmf U 
