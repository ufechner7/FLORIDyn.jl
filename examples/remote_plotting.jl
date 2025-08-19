# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

if Threads.nthreads() > 1
    function init_plotting()
        # Only add a worker if we don't have any dedicated worker processes
        if nprocs() < 2  # nprocs() counts main + workers, so < 2 means no dedicated workers
            println("No dedicated workers found, adding 1 worker...")
            addprocs(1)
            if workers() != [2]
                sleep(1)
            end
            @assert workers() == [2]  # Ensure we have exactly one worker now
            @spawnat 2 eval(:(using ControlPlots))
            @eval @everywhere using FLORIDyn      # Ensure FLORIDyn (including WindFarm) is available on all workers
            
            # Create a completely isolated plt instance for this specific task
            @everywhere function rmt_plot_flow_field(wf, X, Y, Z, vis; msr=EffWind, fig=nothing)
                # Create a fresh plt instance just for this task
                local_plt = ControlPlots.plt
                return plotFlowField(local_plt, wf, X, Y, Z, vis; msr, fig)
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
                # Create a fresh plt instance just for this task
                local_plt = ControlPlots.plt
                return plotMeasurements(local_plt, wf, md, vis; separated=separated, msr)
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
