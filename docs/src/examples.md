# Examples

```@meta
CurrentModule = FLORIDyn
```

## Nordsee One

Below is the wind farm layout used in the Nordsee One example:

![Windfarm layout: Nordsee One](windfarm-layout-nordsee-one.webp)

Source: [https://www.nordseeone.com/windfarm](https://www.nordseeone.com/windfarm)

The configuration file is `2021_54T_NordseeOne.yaml` in the data folder.


### Velocity reduction visualization

The following figure shows the velocity reduction field produced by the example workflow:

![Velocity reduction](ff_velocity_reduction.png)

The example assumes an inflow with constant wind speed and no turbulence. The wind direction changes during the simulation from `255°` to `195°`. The duration of the simulation was 20 minutes. The wind direction is the direction where the wind is coming from, clockwise positive, with `0°` defined as wind from the North.


