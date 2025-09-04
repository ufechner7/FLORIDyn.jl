# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

if Threads.nthreads() > 1
    function init_plotting()
        # Only add a worker if we don't have any dedicated worker processes
        if nprocs() < 2  # nprocs() counts main + workers, so < 2 means no dedicated workers
            println("No dedicated workers found, adding 1 worker...")
            if workers() < [2]
                addprocs(1; exeflags=["-t 1", "--project", "--gcthreads=1,0"])
            end
            if workers() != [2]
                sleep(1)
                @warn "workers: $(workers()), nthreads: $(Threads.nthreads())"
            end
            @assert workers() == [2]  # Ensure we have exactly one worker now
            @spawnat 2 eval(:(using ControlPlots))
            @eval @everywhere using FLORIDyn      # Ensure FLORIDyn (including WindFarm) is available on all workers
            
            # Use a completely isolated plt instance for this specific task
            @everywhere function rmt_plot_flow_field(wf, X, Y, Z, vis; msr=EffWind, fig=nothing)
                local_plt = ControlPlots.plt
                plotFlowField(local_plt, wf, X, Y, Z, vis; msr, fig)
                nothing
            end
            @everywhere function rmt_plot_flow_field(wf, X, Y, Z, vis, t_rel; msr=VelReduction)
                global plot_state
                if abs(t_rel) < 1e-6
                    plot_state = nothing
                end
                local_plt = ControlPlots.plt
                plot_state = plotFlowField(plot_state, local_plt, wf, X, Y, Z, vis, t_rel; msr)
                nothing
            end
            @everywhere function rmt_plot_measurements(wf, md, vis; separated, msr=VelReduction)
                # Use a fresh plt instance just for this task
                local_plt = ControlPlots.plt
                # Pass pltctrl=ControlPlots so that any internal plot_x calls on the (single-threaded) worker
                # have the ControlPlots module available for plotting.
                plotMeasurements(local_plt, wf, md, vis; separated=separated, msr, pltctrl=ControlPlots)
                nothing
            end
            @everywhere function rmt_plotx(times, plot_data...; ylabels=nothing, labels=nothing,
                                            fig="Wind Direction", xlabel="rel_time [s]", ysize=10, bottom=0.02, 
                                            legend_size=nothing, loc=nothing)
                p=ControlPlots.plotx(times, plot_data...; ylabels, labels, fig=fig, xlabel, ysize, bottom, 
                                     legend_size, loc)
                display(p)  # Ensure the plot is displayed
                nothing
            end
            @everywhere function rmt_plot(X, Ys; xlabel, ylabel, labels, xlims, ylims, ann, scatter, title, fig, ysize)
                if isnothing(labels)
                    p = pltctrl.plot(X, Ys; xlabel, ylabel, xlims, ylims, ann, scatter, title, fig, ysize)
                else
                    p = pltctrl.plot(X, Ys; xlabel, ylabel, labels, xlims, ylims, ann, scatter, title, fig, ysize)
                end
                display(p)  # Ensure the plot is displayed
                nothing
            end
            @everywhere function rmt_close_all()
                local_plt = ControlPlots.plt
                return local_plt.close("all")
            end         
        else
            println("Already have $(nprocs()-1) dedicated worker(s), total processes: $(nprocs())")
            println("Workers: $(workers())")
            
            # Ensure ControlPlots and FLORIDyn are loaded on existing workers
            @spawnat 2 eval(:(using ControlPlots))
            @eval @everywhere using FLORIDyn
        end
        nothing
    end
    init_plotting()  # This sets up workers and remote plotting capabilities   
end
