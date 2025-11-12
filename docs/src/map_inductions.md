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
At the beginning, the demand is 40%. At $t=t_\text{start}$ the demand starts to rise, and at $t=t_\text{end}$ it reaches 
80% and then stays constant. In the first test case we use $t_\text{start}=1740s$ and $t_\text{end}=2460s$. The resulting
demand profile is shown in Fig. \ref{fig:power-demand-1t}. 

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

This would be correct without wakes. Because of the wake effects, we need to increase the
set-power of the turbines with a correction factor. This - time dependent - correction factor
is defined as a cubic Hermite spline, based on $n$ control points. The first control point defines the correction 
for $t <= t_\text{start}$, the last control point defines the correction for $t >= t_\text{end}$, and the additional 
control points are distributed evenly between the first and the last point.

The correction function is defined as a piecewise cubic Hermite spline:
\begin{equation}
\label{eq:correction-func}
c(t) = \begin{cases}
c_1 & \text{if } t \leq t_{\text{start}} \\
H_i(s) & \text{if } t_{\text{start}} < t < t_{\text{end}} \\
c_n & \text{if } t \geq t_{\text{end}}
\end{cases}
\end{equation}

where $s = \frac{t - t_{\text{start}}}{t_{\text{end}} - t_{\text{start}}}$ is the normalized time parameter $s \in [0, 1]$, and $H_i(s)$ is the cubic Hermite interpolation between control points $c_i$ at positions $s_i = \frac{i-1}{n-1}$ for $i = 1, \ldots, n$.

For each segment between control points, the Hermite spline is given by:
\begin{equation}
\label{eq:hermite-spline}
H_i(s) = h_{00}(u) \cdot c_i + h_{10}(u) \cdot m_i \cdot \Delta s + h_{01}(u) \cdot c_{i+1} + h_{11}(u) \cdot m_{i+1} \cdot \Delta s
\end{equation}

where $u = \frac{s - s_i}{\Delta s}$ with $\Delta s = \frac{1}{n-1}$, and the Hermite basis functions are:
\begin{align}
h_{00}(u) &= 2u^3 - 3u^2 + 1 \\
h_{10}(u) &= u^3 - 2u^2 + u \\
h_{01}(u) &= -2u^3 + 3u^2 \\
h_{11}(u) &= u^3 - u^2
\end{align}

The tangents $m_i$ at each control point are computed using central differences:
\begin{equation}
\label{eq:tangents}
m_i = \begin{cases}
\frac{c_2 - c_1}{\Delta s} & \text{if } i = 1 \\
\frac{c_{i+1} - c_{i-1}}{2\Delta s} & \text{if } 1 < i < n \\
\frac{c_n - c_{n-1}}{\Delta s} & \text{if } i = n
\end{cases}
\end{equation}

If we combine Eq. \ref{eq:induction-time} and Eq. \ref{eq:correction-func}, we get
\begin{equation}
\label{eq:induction-time-corrected}
a = f(c(t) * d(t))
\end{equation}


Fig. \ref{fig:power-demand-1t} shows the resulting induction factor as a function of time:

![Relative Power and Demand\label{fig:power-demand-1t}](Rel_Power_and_Demand_1T.png){width=70%}



\newpage

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