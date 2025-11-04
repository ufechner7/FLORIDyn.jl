# Map free variables to the induction factors

## Collective turbine control
The most simple way to match production and demand is to control the induction factor of the turbines of the wind farm cluster, assuming all use the same induction factor.

If the free-flow wind velocity $v_w$ is known, we can write:
$$
a = f(t)
$$

## Turbine group (TG) control
The axial induction factor $a$ of each turbine group shall be controlled to achieve the best match between power demand and power production.

Let the relative power demand be defined by the function:
$$
d = g(t) \\
0 < d < 1.0
$$

We can write:
$$
a_i = f(i, t) \\

i \in 1..n
$$