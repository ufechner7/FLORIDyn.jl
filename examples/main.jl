## MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
# Link folder with simulation data

using FLORIDyn, YAML

# get the settings for the wind field and the simulator
wind, sim, con = setup("data/2021_9T_Data.yaml")

# create settings struct
set = Settings(wind)

