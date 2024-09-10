using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf


function run()
    @load("data_tpsci_04.jld2")


    display(v0a)
    n1, n2, sz1, sz2 = correlation_functions(v0a)

    nroots = length(n1)
    
    max_val = 0    
    for r in 1:nroots
        max_val = max(max_val, maximum(abs.(n2[r])))
    end

    for r in 1:nroots
        ordering = [1,3,4,5,2]
        n2r = n2[r][ordering,:][:,ordering]
        display(n2r)
        plotd = heatmap(n2r; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,5.5), ylims = (0.5,5.5),
        yflip=true)
        m, n = size(n2r)
        #vline!(0.5:(n+0.5), c=:white, label=false)
        #hline!(0.5:(m+0.5), c=:white, label=false)
        vline!(0.5:(n+0.5), c=:grey, label=false)
        hline!(0.5:(m+0.5), c=:grey, label=false)


        savefig(plotd,@sprintf("n_correlation_%1i.png", r))
    end
    

    max_val = 0    
    for r in 1:nroots
        max_val = max(max_val, maximum(abs.(sz2[r])))
    end

    for r in 1:nroots
        ordering = [1,3,4,5,2]
        n2r = sz2[r][ordering,:][:,ordering]
        display(n2r)
        plotd = heatmap(n2r; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,5.5), ylims = (0.5,5.5),
        yflip=true)
        m, n = size(n2r)
        vline!(0.5:(n+0.5), c=:grey, label=false)
        hline!(0.5:(m+0.5), c=:grey, label=false)


        savefig(plotd,@sprintf("sz_correlation_%1i.png", r))
    end
end

run()


