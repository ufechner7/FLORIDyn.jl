using Test
using FLORIDyn
using Pkg

@testset "Copy Functions Tests" begin
    # Create a temporary directory for testing
    test_dir = mktempdir()
    original_dir = pwd()
    
    try
        cd(test_dir)
        
        @testset "copy_bin function" begin
            # Test that copy_bin creates bin directory and copies run_julia script
            @test !isdir("bin")
            
            # Call copy_bin function (not exported, so call directly)
            FLORIDyn.copy_bin()
            
            # Check that bin directory was created
            @test isdir("bin")
            
            # Check that run_julia script was copied
            @test isfile("bin/run_julia")
            
            # Check that the file has executable permissions
            stat_info = stat("bin/run_julia")
            @test (stat_info.mode & 0o111) != 0  # Check if any execute bits are set
            
            # Check that the file content is copied correctly
            # Read the original file from the package
            pkg_path = pkgdir(FLORIDyn)
            original_file = joinpath(pkg_path, "bin", "run_julia")
            if isfile(original_file)
                original_content = read(original_file, String)
                copied_content = read("bin/run_julia", String)
                @test original_content == copied_content
            end
            
            # Test that calling copy_bin again doesn't fail (should overwrite)
            @test_nowarn FLORIDyn.copy_bin()
            @test isfile("bin/run_julia")
        end
        
        @testset "copy_bin error handling" begin
            # Test behavior when source package path cannot be found
            # This is hard to test directly, but we can at least ensure
            # the function doesn't crash unexpectedly
            @test_nowarn FLORIDyn.copy_bin()
        end
        
        @testset "copy_bin with existing bin directory" begin
            # Test when bin directory already exists
            mkpath("bin2")
            cd("bin2")
            mkdir("bin")
            
            # Should still work and copy the file
            @test_nowarn FLORIDyn.copy_bin()
            @test isfile("bin/run_julia")
        end
        
    finally
        cd(original_dir)
        # Clean up: remove the temporary directory
        rm(test_dir, recursive=true, force=true)
    end
end

@testset "copy_bin permissions test" begin
    # Test file permissions more thoroughly
    test_dir = mktempdir()
    original_dir = pwd()
    
    try
        cd(test_dir)
        
        FLORIDyn.copy_bin()
        
        if isfile("bin/run_julia")
            stat_info = stat("bin/run_julia")
            # Check that the file has the expected permissions (0o774)
            # On some systems, umask might affect this, so we check for reasonable permissions
            mode = stat_info.mode & 0o777
            
            # Should have read and execute for owner and group, at minimum
            @test (mode & 0o500) == 0o500  # Owner read and execute
            @test (mode & 0o050) == 0o050  # Group read and execute
            
            # Should be executable by owner
            @test (mode & 0o100) != 0
        end
        
    finally
        cd(original_dir)
        rm(test_dir, recursive=true, force=true)
    end
end

@testset "copy_bin integration test" begin
    # Test that the copied script is actually functional
    test_dir = mktempdir()
    original_dir = pwd()
    
    try
        cd(test_dir)
        
        FLORIDyn.copy_bin()
        
        if isfile("bin/run_julia")
            # Test that the script exists and has a shebang or julia reference
            content = read("bin/run_julia", String)
            @test occursin("julia", lowercase(content)) || occursin("#!/", content)
        end
        
    finally
        cd(original_dir)
        rm(test_dir, recursive=true, force=true)
    end
end
