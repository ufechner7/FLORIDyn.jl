# Simple performance comparison between old and new threading implementations
using FLORIDyn
using BenchmarkTools
using OhMyThreads

# Create a simple test case
function test_tmap_performance()
    # Simple computation to test tmap vs traditional approaches
    data = 1:100  # Smaller range to avoid issues
    
    # Using OhMyThreads.tmap with dynamic scheduling
    result_tmap = @btime tmap($data; scheduler=:dynamic) do x
        x^2 + 1  # Simple computation
    end
    
    # Using regular map (serial)
    result_map = @btime map($data) do x
        x^2 + 1  # Simple computation
    end
    
    println("✓ OhMyThreads.tmap result length: ", length(result_tmap))
    println("✓ Regular map result length: ", length(result_map))
    println("✓ Results are equal: ", result_tmap == result_map)
end

test_tmap_performance()
