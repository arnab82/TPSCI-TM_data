using LinearAlgebra, Plots, Statistics

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
J12_=4.776039494
J13_=4.472373054
J14_=4.535551076
J23_=4.532296957
J24_=4.685875743
J34_=4.728174757

H1=zeros(Float64, 3^4, 3^4)
H1=-J12_*kron(Sp_site,Sm_site,I3,I3)-J12_*kron(Sm_site,Sp_site,I3,I3)-J12_*2*kron(Sz_site,Sz_site,I3,I3)-J13_*kron(Sp_site,I3,Sm_site,I3)-J13_*kron(Sm_site,I3,Sp_site,I3)-J13_*2*kron(Sz_site,I3,Sz_site,I3)  -J14_*kron(Sp_site,I3,I3,Sm_site)-J14_*kron(Sm_site,I3,I3,Sp_site)-J14_*2*kron(Sz_site,I3,I3,Sz_site)-J23_*kron(I3,Sp_site,Sm_site,I3)-J23_*kron(I3,Sm_site,Sp_site,I3)-J23_*2*kron(I3,Sz_site,Sz_site,I3)-J24_*kron(I3,Sp_site,I3,Sm_site)-J24_*kron(I3,Sm_site,I3,Sp_site)-J24_*2*kron(I3,Sz_site,I3,Sz_site)-J34_*kron(I3,I3,Sp_site,Sm_site)-J34_*kron(I3,I3,Sm_site,Sp_site)-J34_*2*kron(I3,I3,Sz_site,Sz_site)
H2=zeros(Float64, 3^4, 3^4)
H2=-J12*kron(Sp_site,Sm_site,I3,I3)-J12*kron(Sm_site,Sp_site,I3,I3)-J12*2*kron(Sz_site,Sz_site,I3,I3)-J13*kron(Sp_site,I3,Sm_site,I3)-J13*kron(Sm_site,I3,Sp_site,I3)-J13*2*kron(Sz_site,I3,Sz_site,I3)  -J14*kron(Sp_site,I3,I3,Sm_site)-J14*kron(Sm_site,I3,I3,Sp_site)-J14*2*kron(Sz_site,I3,I3,Sz_site)-J23*kron(I3,Sp_site,Sm_site,I3)-J23*kron(I3,Sm_site,Sp_site,I3)-J23*2*kron(I3,Sz_site,Sz_site,I3)-J24*kron(I3,Sp_site,I3,Sm_site)-J24*kron(I3,Sm_site,I3,Sp_site)-J24*2*kron(I3,Sz_site,I3,Sz_site)-J34*kron(I3,I3,Sp_site,Sm_site)-J34*kron(I3,I3,Sm_site,Sp_site)-J34*2*kron(I3,I3,Sz_site,Sz_site)

# Function to group eigenvalues
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

# Compute eigenvalues for Variational and Extrapolated cases
F_var = eigen(H1)  # Variational
F_ext = eigen(H2)  # Extrapolated

eigenvalues_var = sort(real.(F_var.values))
eigenvalues_ext = sort(real.(F_ext.values))

# Shifted eigenvalues for better visualization
shifted_eigenvalues_var = eigenvalues_var .- minimum(eigenvalues_var)
shifted_eigenvalues_ext = eigenvalues_ext .- minimum(eigenvalues_ext)

# Grouping eigenvalues
grouped_var = group_eigenvalues(shifted_eigenvalues_var)
grouped_ext = group_eigenvalues(shifted_eigenvalues_ext)

# Extract mean, std, and counts for both cases
y_var = [mean(group) for group in grouped_var]
yerr_var = [std(group) for group in grouped_var]
counts_var = [length(group) for group in grouped_var]

y_ext = [mean(group) for group in grouped_ext]
yerr_ext = [std(group) for group in grouped_ext]
counts_ext = [length(group) for group in grouped_ext]

# Define common x-axis (S² values)
xticks = [20,12,12,12,6,6,6,6,6,6,2,2,2,2,2,2,0,0,0]
x_base = 1:length(xticks)

# Adjust x-axis for side-by-side bars
bar_width = 0.3  # Adjust width
x_var = x_base .- bar_width  # Shift left for Variational
x_ext = x_base .+ bar_width  # Shift right for Extrapolated

# Create side-by-side bar plots (NOT overlaying)
p = bar(x_var, y_var, yerr=yerr_var, label="Variational", color=:blue, 
        alpha=0.8, legend=:topleft, width=bar_width, xticks=(x_base, string.(xticks)))

bar!(x_ext, y_ext, yerr=yerr_ext, label="Extrapolated", color=:red, 
     alpha=0.8, width=bar_width)
# Labels
xlabel!("S² value")
ylabel!("ΔE (wavenumber)")
# title!("Comparison of Variational and Extrapolated Energies")

# Save and display
display(p)
savefig(p, "eigenvalues_comparison.pdf")
png(p, "eigenvalues_comparison.png")