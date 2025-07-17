## MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
# Link folder with simulation data

using FLORIDyn, YAML

YAML.load_file("data/2021_9T_Data.yaml")

# Get the settings for the wind field, visualization, controller and Sim.
# Wind, Vis, Sim, Con = setup("2021_9T_Data")
