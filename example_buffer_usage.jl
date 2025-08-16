# Example demonstrating the buffer struct usage for runFLORIS

using FLORIDyn

# Example of how to use the new RunFLORISBuffers for improved performance
function example_usage()
    println("=== RunFLORISBuffers Example ===")
    
    # Create a buffer for 10 rotor points (adjust based on your discretization)
    n_points = 10
    buffers = FLORIDyn.RunFLORISBuffers(n_points)
    
    println("Created buffers for $n_points rotor points:")
    println("  tmp_RPs size: $(size(buffers.tmp_RPs))")
    println("  Vector buffer length: $(length(buffers.cw_y))")
    println("  Buffer types:")
    println("    Float64 arrays: tmp_RPs, cw_y, cw_z, phi_cw, r_cw, tmp_RPs_r, gaussAbs, gaussWght, exp_y, exp_z")
    println("    Bool arrays: core, nw, fw, not_core")
    
    println("\nFor optimal performance:")
    println("1. Create buffers once: buffers = RunFLORISBuffers(n_points)")
    println("2. Reuse in multiple calls: runFLORIS(buffers, set, location_t, states_wf, states_t, d_rotor, floris, windshear)")
    println("3. Or use backward-compatible wrapper: runFLORIS(set, location_t, states_wf, states_t, d_rotor, floris, windshear)")
    
    return buffers
end

if abspath(PROGRAM_FILE) == @__FILE__
    example_usage()
end
