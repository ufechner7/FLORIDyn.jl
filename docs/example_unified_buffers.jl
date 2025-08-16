# Example demonstrating the use of UnifiedBuffers for memory-efficient flow field calculations

using FLORIDyn

# Set up a simple wind farm using default settings
settings_file = "data/2021_9T_Data.yaml"
wind, sim, con, floris, floridyn, ta = setup(settings_file)
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Create a single unified buffer that contains all arrays needed by both
# interpolateOPs! and setUpTmpWFAndRun! functions
unified_buffers = create_unified_buffers(wf)

# The unified_buffers struct contains all the pre-allocated arrays:
# - dist_buffer: for distance calculations in interpolateOPs!
# - sorted_indices_buffer: for sorting in interpolateOPs!
# - M_buffer: main result buffer for setUpTmpWFAndRun!
# - iTWFState_buffer: turbine wind field state buffer
# - tmp_Tpos_buffer: temporary turbine position buffer
# - tmp_WF_buffer: temporary wind field buffer
# - tmp_Tst_buffer: temporary turbine state buffer
# - dists_buffer: distance buffer for setUpTmpWFAndRun!
# - plot_WF_buffer: wind field plotting buffer
# - plot_OP_buffer: observation point plotting buffer

println("UnifiedBuffers created successfully!")
println("Type: ", typeof(unified_buffers))
println("Available buffer fields:")
for field in fieldnames(typeof(unified_buffers))
    buffer = getfield(unified_buffers, field)
    println("  $field: $(typeof(buffer)) of size $(size(buffer))")
end

# This approach reduces memory allocations by:
# 1. Creating a single struct containing all buffers
# 2. Reusing the same buffers across multiple function calls
# 3. Eliminating the need to pass individual buffer arrays separately
# 4. Providing a clean, type-safe interface for buffer management

# The getMeasurementsP function now uses this unified approach internally,
# providing better performance with reduced memory allocations.
