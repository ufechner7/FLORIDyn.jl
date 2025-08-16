using FLORIDyn
using Test

@testset "UnifiedBuffers Integration Tests" begin
    
    @testset "create_unified_buffers basic functionality" begin
        # Create a minimal wind farm for testing
        wf = WindFarm(
            nT = 2,
            nOP = 10,
            posBase = [0.0 1000.0; 0.0 0.0],
            D = [120.0, 120.0],
            Names_T = ["T1", "T2"],
            Names_WF = ["U_eff", "TI_eff"],
            Names_OP = ["U", "Dir", "TI"],
            States_WF = rand(5, 2),  # 5 timesteps × 2 WF variables
            States_T = rand(3, 2),   # 3 states × 2 turbines
            States_OP = rand(10, 3), # 10 observation points × 3 OP variables
            dep = [Int[], Int[]]     # No dependencies initially
        )
        
        # Test default method
        buffers = create_unified_buffers(wf)
        @test typeof(buffers) == UnifiedBuffers
        @test !isnothing(buffers.floris_buffers)
        
        # Test with rotor points parameter  
        buffers_with_points = create_unified_buffers(wf, 5)
        @test typeof(buffers_with_points) == UnifiedBuffers
        @test !isnothing(buffers_with_points.floris_buffers)
        
        # Verify buffer sizes
        @test length(buffers.dist_buffer) == wf.nOP
        @test length(buffers.sorted_indices_buffer) == wf.nOP
        @test size(buffers.tmp_Tpos_buffer, 1) == wf.nT + 1
        @test size(buffers.tmp_WF_buffer, 1) == wf.nT + 1
        @test size(buffers.tmp_Tst_buffer, 1) == wf.nT + 1
        @test length(buffers.dists_buffer) == wf.nT + 1
    end
    
    @testset "FLORIS buffers integration" begin
        # Create wind farm and FLORIS parameters
        wf = WindFarm(
            nT = 2,
            nOP = 10,
            posBase = [0.0 1000.0; 0.0 0.0],
            D = [120.0, 120.0],
            Names_T = ["T1", "T2"],
            Names_WF = ["U_eff", "TI_eff"],
            Names_OP = ["U", "Dir", "TI"],
            States_WF = rand(5, 2),
            States_T = rand(3, 2),
            States_OP = rand(10, 3),
            dep = [Int[], Int[]]
        )
        floris = Floris(
            alpha = 2.32,
            beta = 0.154,
            k_a = 0.3837,
            k_b = 0.0037,
            k_fa = 0.73,
            k_fb = 0.8325,
            k_fc = 0.0325,
            k_fd = -0.32,
            eta = 1,
            p_p = 2.2,
            airDen = 1.225,
            TIexp = 3,
            rotor_points = 1
        )
        
        # Test FLORIS method dispatch
        buffers = create_unified_buffers(wf, floris)
        @test typeof(buffers) == UnifiedBuffers
        @test !isnothing(buffers.floris_buffers)
        @test typeof(buffers.floris_buffers) == FLORIDyn.RunFLORISBuffers
        
        # Verify FLORIS buffer fields exist
        fb = buffers.floris_buffers
        @test hasfield(typeof(fb), :tmp_RPs)
        @test hasfield(typeof(fb), :cw_y)
        @test hasfield(typeof(fb), :cw_z)
        @test hasfield(typeof(fb), :phi_cw)
        @test hasfield(typeof(fb), :r_cw)
        @test hasfield(typeof(fb), :core)
        @test hasfield(typeof(fb), :nw)
        @test hasfield(typeof(fb), :fw)
        @test hasfield(typeof(fb), :tmp_RPs_r)
        @test hasfield(typeof(fb), :gaussAbs)
        @test hasfield(typeof(fb), :gaussWght)
        @test hasfield(typeof(fb), :exp_y)
        @test hasfield(typeof(fb), :exp_z)
        @test hasfield(typeof(fb), :not_core)
    end
    
    @testset "Buffer allocation check" begin
        # Test that buffers can handle allocation tracking
        wf = WindFarm(
            nT = 2,
            nOP = 10,
            posBase = [0.0 1000.0; 0.0 0.0],
            D = [120.0, 120.0],
            Names_T = ["T1", "T2"],
            Names_WF = ["U_eff", "TI_eff"],
            Names_OP = ["U", "Dir", "TI"],
            States_WF = rand(5, 2),
            States_T = rand(3, 2),
            States_OP = rand(10, 3),
            dep = [Int[], Int[]]
        )
        
        buffers = create_unified_buffers(wf, 10)
        fb = buffers.floris_buffers
        
        # Check that all arrays are properly allocated
        @test length(fb.tmp_RPs) > 0
        @test length(fb.cw_y) > 0
        @test length(fb.cw_z) > 0
        @test length(fb.phi_cw) > 0
        @test length(fb.r_cw) > 0
        @test length(fb.core) > 0
        @test length(fb.nw) > 0
        @test length(fb.fw) > 0
        @test length(fb.tmp_RPs_r) > 0
        @test length(fb.gaussAbs) > 0
        @test length(fb.gaussWght) > 0
        @test length(fb.exp_y) > 0
        @test length(fb.exp_z) > 0
        @test length(fb.not_core) > 0
    end
end
