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

