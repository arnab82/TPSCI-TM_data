using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf
using ActiveSpaceSolvers

C = npzread("mo_coeffs.npy")
h0 = npzread("ints_h0.npy")
h1 = npzread("ints_h1.npy")
h2 = npzread("ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

Pa = npzread("Pa.npy")
Pb = npzread("Pb.npy")
@printf(" Input energy:    %12.8f\n", compute_energy(ints, RDM1(Pa, Pb)))


init_fspace =  [(5, 4), (5, 0), (3, 3)]
clusters    =  [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]]
clusters = [MOCluster(i, collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

rdm1 = RDM1(n_orb(ints))


# # Do CMF
ansatze = [FCIAnsatz(10, 5,4), FCIAnsatz(10,5,0),FCIAnsatz(12,3,3)]
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo_newton(ints, clusters, init_fspace,ansatze,rdm1, maxiter_oo = 400,
                           tol_oo=1e-6, 
                           tol_d1=1e-9, 
                           tol_ci=1e-11,
                           step_trust_region=1.5,
                           verbose=0, 
                            trust_region=true,
                           sequential=true)
ints = orbital_rotation(ints, U)
C= C*U
npzwrite("Ccmf_32.npy", C)
@save "data_cmf_32_Mn_Cy.jld2" clusters init_fspace ints d1 e_cmf U 
