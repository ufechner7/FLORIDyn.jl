using DistributedNext

# Only add a worker if we don't have any workers yet
if nworkers() == 0
    println("No workers found, adding 1 worker...")
    addprocs(1)
else
    println("Already have $(nworkers()) worker(s), using existing workers")
end

# Define function on all processes (main + workers)
@everywhere function identity_func(x)
    return x
end

function pltr(data)
    worker_id = workers()[1]
    remote_do(identity_func, worker_id, data)
end

data = rand(UInt8, 1024^2)  # 1MB random data
worker_id = workers()[1]  # Get the first available worker ID
println("Using worker $worker_id")

# Pre-compile by calling once
pltr(data)
@time pltr(data)
nothing
