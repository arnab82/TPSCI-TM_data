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

    e = e .- (sum(e)/length(e))

    S = vproj'*vproj
    println(" S:")
    display(S)
    vproj_orth = vproj * inv(sqrt(S))

    Heff = vproj_orth * diagm(e) * vproj_orth'

    println(" Heff")
    display(Heff)

    # shift diagonal to center at zero
    return Heff
end

@load("data_tpsci_04.jld2")

conversion = 219474.63

display(FermiCG.get_vector(v0))
v_bare = v0
e_bare = eci
e_exact = (e2a + e0a)
v_exact = deepcopy(v0a)
v_model = deepcopy(v0)
v_diabat = deepcopy(v_model)
FermiCG.eye!(v_diabat)

# First make bare 
println()
println(" Make Model Space Bare Hamiltonian")
Heff_bare = make_H(v_diabat, v_bare, e_bare)
for i in 1:length(e_bare)
    @printf(" Adiabatic: %12.8f Diabatic: %12.8f\n", e_bare[i], Heff_bare[i,i])
end
Heff_bare .*= conversion

# Now exact
println()
println(" Make Model Space Effective Hamiltonian")
Heff_exact = make_H(v_diabat, v_exact, e_exact)
for i in 1:length(e_exact)
    @printf(" Adiabatic: %12.8f Diabatic: %12.8f\n", e_exact[i], Heff_exact[i,i])
end
Heff_exact .*= conversion

max_val = max(maximum(abs.(Heff_bare)), maximum(abs.(Heff_exact)))

#plotd = heatmap(Heff*conversion; color=palette([:teal, :white, :orange], 10), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm)
plotd = heatmap(Heff_bare; color=palette(:RdGy_10, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                xlims = (0.5,4.5), ylims = (0.5,4.5),
                yflip=true)
m, n = size(Heff_bare)
vline!(0.5:(n+0.5), c=:white, label=false)
hline!(0.5:(m+0.5), c=:white, label=false)


savefig(plotd,"Heff_bare.png")



plotd = heatmap(Heff_exact; color=palette(:RdGy_10, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                clims=(-max_val, max_val), ticks = false, xaxis=false,yaxis=false,
                xlims = (0.5,4.5), ylims = (0.5,4.5),
                yflip=true)
#plotd = heatmap(Heff*27.2114079527*1000; color=palette(:Blues, 10), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm, xlims = (0,5), ylims = (0,5))
m, n = size(Heff_bare)
vline!(0.5:(n+0.5), c=:white, label=false)
hline!(0.5:(m+0.5), c=:white, label=false)
savefig(plotd,"Heff_exact.png")


v_bare_mat = FermiCG.get_vector(v0)
S2 = v_bare_mat * diagm([12, 6, 2, 0] ) * v_bare_mat'
println(" S2: ")
display(S2)
