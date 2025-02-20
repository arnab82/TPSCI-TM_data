using LinearAlgebra


I3 = Matrix{Float64}(I, 3, 3)
Sz_site = [1 0 0; 0 0 0; 0 0 -1]
Sp_site = √2 * [0 1 0; 0 0 1; 0 0 0]
Sm_site = √2 * [0 0 0; 1 0 0; 0 1 0]
# EXTRAPOLATION
J12=	4.633313611
J13=	4.237422671
J14=	4.407452587
J23=	4.363819426
J24=	4.47820604
J34=	4.606011794
#VARIATIONAL
# J12=4.776039494
# J13=4.472373054
# J14=4.535551076
# J23=4.532296957
# J24=4.685875743
# J34=4.728174757

H=zeros(Float64, 3^4, 3^4)
H=-J12*kron(Sp_site,Sm_site,I3,I3)-J12*kron(Sm_site,Sp_site,I3,I3)-J12*2*kron(Sz_site,Sz_site,I3,I3)-J13*kron(Sp_site,I3,Sm_site,I3)-J13*kron(Sm_site,I3,Sp_site,I3)-J13*2*kron(Sz_site,I3,Sz_site,I3)  -J14*kron(Sp_site,I3,I3,Sm_site)-J14*kron(Sm_site,I3,I3,Sp_site)-J14*2*kron(Sz_site,I3,I3,Sz_site)-J23*kron(I3,Sp_site,Sm_site,I3)-J23*kron(I3,Sm_site,Sp_site,I3)-J23*2*kron(I3,Sz_site,Sz_site,I3)-J24*kron(I3,Sp_site,I3,Sm_site)-J24*kron(I3,Sm_site,I3,Sp_site)-J24*2*kron(I3,Sz_site,I3,Sz_site)-J34*kron(I3,I3,Sp_site,Sm_site)-J34*kron(I3,I3,Sm_site,Sp_site)-J34*2*kron(I3,I3,Sz_site,Sz_site)


# Calculate and print the eigenvalues
F = eigen(H)
eigenvalues = F.values
F.vectors
eigenvalues = eigvals(H)

display(H)
sorted_eigenvalues = sort(real.(eigenvalues))
shifted_eigenvalues = sorted_eigenvalues .- minimum(sorted_eigenvalues)
println(sort(real.(eigenvalues)))
println("Lowest 19 eigenvalues:")
display(sort(real.(eigenvalues)))

# Calculate the ground state

ground_state = eigvecs(Hermitian(H))[:, 1]
# for i in 1:3^4
#     println("state component $i: ", state[i])
#     state[i] = eigvecs(Hermitian(H))[:, i]
# end

using Plots
using Statistics


function group_eigenvalues(eigenvalues, tolerance=1e-6)
    grouped = []
    current_group = [eigenvalues[1]]
    
    for i in 2:length(eigenvalues)
        if isapprox(eigenvalues[i], current_group[1], atol=tolerance)
            push!(current_group, eigenvalues[i])
        else
            push!(grouped, current_group)
            current_group = [eigenvalues[i]]
        end
    end
    
    push!(grouped, current_group)
    return grouped
end
shifted_eigenvalues = group_eigenvalues(shifted_eigenvalues)
# Prepare data for plotting
x = 1:length(shifted_eigenvalues)
y = [mean(group) for group in shifted_eigenvalues]
yerr = [std(group) for group in shifted_eigenvalues]
counts = [length(group) for group in shifted_eigenvalues]

# Define colors for different groups
colors = [:red, fill(:orange, 3)..., fill(:green, 6)..., fill(:blue, 6)..., fill(:purple, 3)...]

# Create the bar plot
xticks=[20,12,12,12,6,6,6,6,6,6,2,2,2,2,2,2,0,0,0]
x=1:length(xticks)

p = bar(x, y, yerr=yerr, label="", title="",
        xlabel="S² value", ylabel="ΔE (wavenumber)",
        color=colors, xticks=(x, string.(xticks)),
        legend=false,size = (600, 400))



# Display the plot
display(p)

# savefig(p, "eigenvalues_barplot_var.pdf")
# png(p, "eigenvalues_barplot_var.png")
savefig(p, "eigenvalues_barplot_EXTRAP.pdf")
png(p, "eigenvalues_barplot_EXTRAP.png")