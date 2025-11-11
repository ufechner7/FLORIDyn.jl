# Map free variables to the induction factors

## Collective turbine control
The most simple way to match production and demand is to control the induction factor of the turbines of the wind farm cluster, assuming all use the same induction factor.

We assume
\begin{equation}
\label{eq:cp-max}
C_{\text{p,max}} = \frac{16}{27}
\end{equation}

Furthermore, we assume that the relative demand (relative to the power at free-flow velocity) is given as
\begin{equation}
\label{eq:demand-rel}
p_{\text{d,rel}} = d(t)
\end{equation}

The power coefficient of an ideal wind turbine is related to the induction factor following
\begin{equation}
\label{eq:cp-induction}
C_{\text{p}} = 4a(1-a)^2
\end{equation}
The inverse relationship $a(C_{\text{p}})$ is obtained by numerically solving this equation for $a \in [0, 1/3]$.

By combining Eqs. \eqref{eq:cp-max}, \eqref{eq:demand-rel} and \eqref{eq:cp-induction}, and applying $C_{\text{p}} = d(t) \cdot C_{\text{p,max}}$ we can write:
\begin{equation}
\label{eq:induction-time}
a = f(d(t))
\end{equation}

Fig. \ref{fig:power-demand-1t} shows the resulting induction factor as a function of time:

![Relative Power and Demand\label{fig:power-demand-1t}](Rel_Power_and_Demand_1T.png)

## Turbine group (TG) control
The axial induction factor $a$ of each turbine group shall be controlled to achieve the best match between power demand and power production.

Let the relative power demand be defined by the function:
\begin{equation}
\label{eq:demand-func}
d = g(t), \quad 0 < d < 1.0
\end{equation}

We can write:
\begin{equation}
\label{eq:induction-group}
a_i = f(i, t), \quad i \in \{1, \ldots, n\}
\end{equation}