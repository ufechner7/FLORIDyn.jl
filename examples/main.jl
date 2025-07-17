## MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
# Link folder with simulation data

using FLORIDyn, YAML

# get the settings for the wind field and the simulator
wind, sim = setup("data/2021_9T_Data.yaml")

# get the settings for the wind field, visualization, controller and Sim.
# Wind, Vis, Sim, Con = setup("2021_9T_Data")
