function init_plotting()
    # Only add a worker if we don't have any dedicated worker processes
    if nprocs() < 2  # nprocs() counts main + workers, so < 2 means no dedicated workers
        println("No dedicated workers found, adding 1 worker...")
        @time addprocs(1)
        @eval @everywhere using ControlPlots  # Ensure ControlPlots is available on all workers
        @eval @everywhere using FLORIDyn      # Ensure FLORIDyn (including WindFarm) is available on all workers
    else
        println("Already have $(nprocs()-1) dedicated worker(s), total processes: $(nprocs())")
        println("Workers: $(workers())")
    end
end