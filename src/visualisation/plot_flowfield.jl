"""
    getMeasurements(X, Y, nM, zh, T, paramFLORIS, Wind, Vis)

Wrapper function to disguise the grid points as turbines with one rotor
point and experience almost the same calculations as the rotor points in
the simulation.

Single thread version

# Arguments
- `X::Matrix`: X-coordinates of grid points
- `Y::Matrix`: Y-coordinates of grid points  
- `nM::Int`: Number of measurements
- `zh::Real`: Hub height
- `T::WindFarm`: Wind farm object containing turbine data
- `paramFLORIS::Floris`: FLORIS model parameters
- `Wind::Wind`: Wind field configuration
- `Vis`: Visualization parameters (currently unused)

# Returns
- `Z::Array{Float64,3}`: 3D array of measurements with dimensions (size(X,1), size(X,2), nM)
"""
function getMeasurements(X, Y, nM, zh, T, paramFLORIS, Wind, Vis)
    sizeX = size(X)
    Z = zeros(sizeX[1], sizeX[2], nM)
    
    # Single-threaded loop (can be parallelized with @threads or Distributed.@distributed)
    for iGP in 1:length(X)
        xGP = X[iGP]
        yGP = Y[iGP]
        
        GPdep = [1:T.nT]

        # Create T equivalent to get interpolated OPs
        # Initialize GP as a similar structure to WindFarm
        GP = deepcopy(T)  # Start with a copy to get the right structure
        
        # Modify the necessary fields
        GP.dep = Vector{Vector{Int64}}(undef, 1)
        GP.dep[1] = GPdep
        GP.StartI = T.StartI
        GP.nT = 1
        GP.States_OP = T.States_OP
        GP.posBase = reshape([xGP, yGP, 0.0], 1, 3)  # Make it a 1×3 matrix
        GP.nOP = T.nOP
        GP.intOPs = interpolateOPs(GP)
        
        GP.posNac = reshape([0.0, 0.0, zh], 1, 3)  # Make it a 1×3 matrix
        GP.States_WF = T.States_WF
        GP.States_T = T.States_T
        GP.D = vcat(T.D, 0.0)  # Append 0 to the diameter array
        
        tmpM, _ = setUpTmpWFAndRun(paramFLORIS, GP, paramFLORIS, Wind)
        
        # Convert linear index to subscripts
        rw, cl = divrem(iGP - 1, sizeX[1])
        rw += 1
        cl += 1
        Z[rw, cl, 1:3] = tmpM
    end
    
    return Z
end