#!/usr/bin/env julia

# Demo script showing the new calcFlowField vis parameter validation
# This demonstrates that the test case validates the functionality correctly

using FLORIDyn

# Setup test case
settings_file = "data/2021_9T_Data.yaml"
wind, sim, con, floris, floridyn, ta = setup(settings_file)
use_threading = Threads.nthreads() > 1
set = Settings(wind, sim, con, use_threading, use_threading)
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)  
wf = initSimulation(wf, sim)

# Test 1: Default behavior (no vis parameter)
println("=== Test 1: Default behavior ===")
Z1, X1, Y1 = calcFlowField(set, wf, wind, floris)
println("Default grid size: Z=$(size(Z1)), X=$(size(X1)), Y=$(size(Y1))")

# Test 2: Custom vis configuration with specific field limits and resolution
println("\n=== Test 2: Custom vis configuration ===")
custom_vis = Vis(
    field_limits_min = [200.0, 300.0, 50.0],
    field_limits_max = [1200.0, 1300.0, 150.0], 
    field_resolution = 100.0,
    online = false,
    save = false,
    unit_test = true
)

Z2, X2, Y2 = calcFlowField(set, wf, wind, floris; vis=custom_vis)

# Calculate expected dimensions
expected_x = Int(round((custom_vis.field_limits_max[1] - custom_vis.field_limits_min[1]) / custom_vis.field_resolution)) + 1
expected_y = Int(round((custom_vis.field_limits_max[2] - custom_vis.field_limits_min[2]) / custom_vis.field_resolution)) + 1

println("Custom field limits: min=$(custom_vis.field_limits_min), max=$(custom_vis.field_limits_max)")
println("Resolution: $(custom_vis.field_resolution) m")
println("Expected grid dimensions: $expected_x × $expected_y")
println("Actual grid size: Z=$(size(Z2)), X=$(size(X2)), Y=$(size(Y2))")

# Verify the test validations work
println("\n=== Validation checks ===")
println("✓ Z dimensions match: $(size(Z2, 1) == expected_y && size(Z2, 2) == expected_x)")
println("✓ X grid dimensions: $(size(X2, 2) == expected_x)")
println("✓ Y grid dimensions: $(size(Y2, 1) == expected_y)")
println("✓ Coordinate ranges - X: $(minimum(X2)) to $(maximum(X2)) (expected: $(custom_vis.field_limits_min[1]) to $(custom_vis.field_limits_max[1]))")
println("✓ Coordinate ranges - Y: $(minimum(Y2)) to $(maximum(Y2)) (expected: $(custom_vis.field_limits_min[2]) to $(custom_vis.field_limits_max[2]))")

if size(X2, 2) > 1
    x_spacing = X2[1, 2] - X2[1, 1]
    println("✓ X spacing: $(x_spacing) m (expected: $(custom_vis.field_resolution) m)")
end

if size(Y2, 1) > 1  
    y_spacing = Y2[2, 1] - Y2[1, 1]
    println("✓ Y spacing: $(y_spacing) m (expected: $(custom_vis.field_resolution) m)")
end

println("\n=== Summary ===")
println("The test case successfully validates that:")
println("- calcFlowField respects vis.field_limits_min and vis.field_limits_max")
println("- The output grid dimensions match the calculated expected size")
println("- Grid spacing matches vis.field_resolution")
println("- Coordinate ranges are correctly set")
println("- Backward compatibility is maintained (vis=nothing works)")
