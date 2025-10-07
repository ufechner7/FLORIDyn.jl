# FLORIDyn v0.5.6 - 2025-07-10
# Added
- example calc_induction.jl
- example step_response.jl
- example_dev mwe_06.jl (developer example)
- example_dev mwe_07.jl (developer example)
- add turbine groups to the yaml file and functions and structs to work with them
- a new control mode can process time depended induction factors for turbine groups
- all correction functions from the Matlab code have been implemented
- the function `cp_fun` has been added that calculates the power coefficient as function of angle of attack and tip speed ratio
- example_dev cp_
- function plot_rmt supports now all plot variants of `ControlPlots.plot`

# FLORIDyn v0.5.5 - 2025-25-09
## Fixed
- Fix .zenodo.json

# FLORIDyn v0.5.4 - 2025-25-09
## Fixed
- Fix .zenodo.json

# FLORIDyn v0.5.3 - 2025-25-09
## Added
- Add power plots to the package itself, function calc_rel_power()
- Add function plot_rmt
## Fixed
- Fix reading Control_YawConstant.csv
## Changed
- Update main_power_plot.jl
- Delete main_power_plot_II.jl

# FLORIDyn v0.5.2
## Added
- Add function getYaw(::Yaw_Constant, ...)
- Increase size of `SOWFA_nacelleYaw.csv` to allow longer simulations (400s longer)
- Add main_power_plot.jl
- Add main_power_plot_II.jl
- Add parameter `simple=false` to function `prepare_large_plot_inputs()`

## Fixed
- Fix y-label of velocity reduction plot
- Fix example main_large_msr.jl

