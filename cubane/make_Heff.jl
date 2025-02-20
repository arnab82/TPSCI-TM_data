using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf

function make_H(diabats::TPSCIstate, adiabats::TPSCIstate, ein)
    vproj = FermiCG.overlap(diabats, adiabats)
    println(" Vproj:")
    display(vproj)

    e = ein

    #e = e .- (sum(e)/length(e))

    S = vproj'*vproj
    println(" S:")
    display(S)
    vproj_orth = vproj * inv(sqrt(S))

    Heff = vproj_orth * diagm(e) * vproj_orth'

    println(" Heff")
    display(Heff)
    show(IOContext(stdout, :compact => false, :precision => 10), "text/plain", Heff)

    # shift diagonal to center at zero
    return Heff
end


#hartree to meV
conversion = 219474.63

# @load "tpsci_out_Ss_1.jld2"
@load "tpsci_out_Ss_1_nonhosvd_2e4.jld2"
@load "guess_tpsci.jld2"
display(FermiCG.get_vector(ci_vector))
# Root       Energy           S2
#      1 -447.67171797  19.99922022
#      2 -447.67155748  12.00050804
#      3 -447.67155730  12.00047443
#      4 -447.67155319  12.00012204
# e0=[-447.67171797, -447.67155748, -447.67155730, -447.67155319]
e0=[-447.6728006,
-447.6726356,
-447.6726314,
-447.6726294]
e0_extrap=[-447.6744277,
-447.6742696,
-447.6742635,
-447.6742629]
ecore=-8604.886065980987

e0_sf=[
-9207.0307957333,
-9207.0307731251,
-9207.0307731227,
-9207.0307411800]
v_barett = v_guess
e_barett = e_guess
e_exacttt = (e0)
# e_exacttt = (e0_extrap)
# e_exacttt = (e0_sf)
#e_exacttt = (e2 + e0)
v_exacttt = deepcopy(v0a)
v_modeltt = deepcopy(v_guess)
v_diabattt = deepcopy(ci_vector)
FermiCG.eye!(v_diabattt)

# First make bare 
println()
println(" Make Model Space Bare Hamiltonian")
Heff_bare = make_H(v_diabattt, v_barett, e_barett)
Heff_bare .*= conversion

# Now exact
println()
println(" Make Model Space Effective Hamiltonian")
Heff_exact = make_H(v_diabattt, v_exacttt, e_exacttt)
Heff_exact .*= conversion


Heff_bare = Heff_bare - diagm(diag(Heff_bare))
Heff_exact = Heff_exact - diagm(diag(Heff_exact))

display(Heff_bare)
display(Heff_exact)
max_val = max(maximum(abs.(Heff_bare)), maximum(abs.(Heff_exact)))

#plotd = heatmap(Heff_bare_tt; color=palette(:bone_1, rev=true), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
plotd = heatmap(Heff_bare; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                xlims = (0.5,4.5), ylims = (0.5,4.5), yflip=true)
m, n = size(Heff_bare)
vline!(0.5:(n+0.5), c=:white, label=false)
hline!(0.5:(m+0.5), c=:white, label=false)
tmp = [0.5, 1.5, 4.5]
vline!(tmp, c=:grey, label=false)
hline!(tmp, c=:grey, label=false)


savefig(plotd,"Heff_bare_spin.png")

#plotd = heatmap(Heff_exact_tt; color=palette(:bone_1, rev=true), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
plotd = heatmap(Heff_exact; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                clims=(-max_val, max_val), ticks = false, xaxis=false,yaxis=false,
                xlims = (0.5,4.5), ylims = (0.5,4.5), yflip=true)
m, n = size(Heff_bare)
vline!(0.5:(n+0.5), c=:white, label=false)
hline!(0.5:(m+0.5), c=:white, label=false)

tmp = [0.5, 1.5, 4.5]
vline!(tmp, c=:grey, label=false)
hline!(tmp, c=:grey, label=false)
savefig(plotd,"Heff_exact_spin.png")



