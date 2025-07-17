## MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
# Link folder with simulation data

using FLORIDyn, YAML

# get the settings for the wind field and the simulator
wind, sim, con = setup("data/2021_9T_Data.yaml")

# % Add according functions to the search path
# addFLORISPaths
    vel_mode = str2type("Velocity_" * wind.input_vel)
    dir_mode = str2type("Direction_" * wind.input_dir)
    turb_mode = str2type("TI_" * wind.input_ti)
    shear_mode = str2type("Shear_" * wind.input_shear)
    # addpath(['Correction/Direction_' Wind.Correction.Dir]);
    # addpath(['Correction/Velocity_' Wind.Correction.Vel]);
    # addpath(['Correction/TI_' Wind.Correction.TI]);

