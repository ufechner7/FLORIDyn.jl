# Map free variables to the induction factors

## Collective turbine control
The most simple way to match production and demand is to control the induction factor of the turbines of the wind farm cluster, assuming all use the same induction factor.

We assume
$$
C_{\text{p,max}} = \frac{16}{27}
$$

Furthermore, we assume that the relative demand (relative to the power at free-flow velocity) is given as
$$
p_{\text{d,rel}} = d(t)
$$

The power coefficient of an ideal wind turbine is related to the induction factor following
$$
C_{\text{p}} = 4a(1-a)^2
$$
The inverse relationship $a(C_{\text{p}})$ is obtained by numerically solving this equation for $a \in [0, 1/3]$.

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