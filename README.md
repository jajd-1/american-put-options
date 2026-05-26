# American put option pricing as a parabolic free boundary problem 

The goal of this project is to formulate the American put option pricing problem mathematically as a parabolic obstacle/free boundary problem, solve the problem numerically and explore how the free boundary is affected by varying the parameters. 

More precisely, we first construct and validate our numerical solution produced by the finite difference method (FDM) against the closed-form solution to the Black-Scholes equation in the case of European put options. We then incorporate the obstacle imposed in the American case, re-solve the problem numerically and extract the free boundary, which determines when it is optimal for the holder to exercise the option. In the process we address numerical challenges occuring near the free boundary. We then study how the free boundary changes as parameters such as the volatility, interest rate and maturity date are varied. 

## European and American put options

### European case

A European put option gives the holder the right (but not the obligation) to sell an asset at a specified strike price $K$ at a specified maturity time $T$. If $S_t$ denotes the price of the asset at time $t \leq T$, then the payoff at maturity is 

$$(K-S_T)^+.$$

The price of the option at times $t < T$ may be considered as a function of $S$ and $t$ and, under certain assumptions, is given by the solution to the Black-Scholes equation 

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS\frac{\partial V}{\partial S} - rV = 0$$

with terminal condition $V(S,T) = (K-S)^+$. Here, $\sigma$ denotes the volatility and $r$ is the risk-free rate. We also have the boundary conditions $\lim_{S\rightarrow \infty} V(S,t) = 0$ (the put option is essentially worthless when the underlying price is far above the strike price) and $V(0,t) = K e^{-r(T-t)}$ (if the stock is worthless then the put option is worth its strike price discounted from maturity back to time $t$). 

Making the substitutions

$$\tau = T - t, \quad u = Se^{r\tau} \quad \text{and}\quad x = \ln\bigg(\frac{S}{K}\bigg) + \bigg(r - \frac{1}{2}\sigma^2\bigg)\tau$$

transforms the Black-Scholes equation into the diffusion equation

$$\frac{\partial u}{\partial \tau} = \frac{1}{2}\sigma^2 \frac{\partial^2 u}{\partial x^2},$$

which can be solved explicitly and transformed back to the original variables, yielding the following price at time $T$ to expiration:

$$K e^{-rT}N(-d_2) - SN(-d_1),$$

where $N$ is the CDF of the standard normal distribution and 

$$d_1 = \frac{\ln\big(\frac{S}{K}\big) + \big(r + \frac{\sigma^2}{2}\big)}{\sigma\sqrt{T}}, \quad d_2 = d_1 - \sigma\sqrt{T}.$$

### American case

An American put option gives the holder the additional flexibility to exercise the option at any time before maturity. The first change required in the PDE formulation is the boundary condition at $S=0$: since one no longer has to wait until maturity, an American put option is simply worth its strike price if the stock price hits zero, i.e. $V(0,t) = K$. The second change is more interesting: rather than just the terminal condition $V(S,T) = (K-S)^+$, the option price must also satisfy

$$V(S,t) \geq (K-S)^+.$$

at all times $t \leq T$. Indeed, if this weren't the case then there would be an arbitrage opportunity: one could buy the asset at price $S$, purchase the put option at price $V$ and immediately exercise it to sell the asset at price $K$, resulting in a profit of $K - S - V > 0$.

Therefore, the price of an American put option is characterised by the constraint $V(S,t) \geq (K-S)^+$ together with the condition that the Black-Scholes equation is satisfied in the region where $V(S,t) > (K-S)^+$, which we refer to as the *continuation region*. In the *exercise region*, where $V(S,t) = (K-S)^+$, it is optimal for the holder to exercise the option. The interface between these two regions is the so-called free boundary, which determines for each time $t$ the threshold asset price below which exercise is optimal and above which the Black-Scholes equation governs the option price. 

We can picture the free boundary as a curve $\{(S,t): S\in[0,\infty), t\in[0,T]\}\subset\mathbb{R}^2$, although stricty speaking this requires one to prove something about the structure and regularity of free boundaries for Lipschitz obstacles. We put aside this analytic aspect of the free boundary problem for now, and consider instead numerical approaches to the problem. 

## Numerical analysis 

### European case

We consider the equation obtained by making the substitution $\tau = T - t$ and considering $V$ as a function of $S$ and $\tau$:

$$\frac{\partial V}{\partial \tau} =  \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS\frac{\partial V}{\partial S}.$$

Note that the terminal condition has now become the initial condition $V(S,0) = (K-S)^+$. The boundary condition at $S=0$ reads $V(0,\tau) = Ke^{-r\tau}$, and to approximate the boundary condition at $S = \infty$, we truncate the $S$ domain to $S\in[0,S_\text{max}]$, where $S_\text{max}$ is some fixed multiple of the strike price at which the option would essentially be worthless in practice (e.g. $S_\text{max}\approx 5K$), and set $V(S_\text{max},\tau) = 0$. Then the parabolic cylinder on which to solve the PDE is just the rectangle $[0, S_\text{max}]\times(0,T]$, which lends itself to the finite difference method. 

We discretise the interval $[0, S_\text{max}]$ into $M+1$ points 

$$S_i = i\Delta S \quad \text{for} \quad i=0,1,\dots,M \quad \text{and} \quad \Delta S = \frac{S_\text{max}}{M}$$

and the interval $[0,T]$ into $N+1$ points 

$$\tau_n = n\Delta \tau \quad \text{for} \quad i=0,1,\dots,N \quad \text{and} \quad \Delta \tau = \frac{T}{N},$$

and define

$$V_i^n = V(S_i, \tau_n).$$

We start with a fully implicit scheme whereby spatial derivatives are evaluated using backward differences at time $t_{n+1}$. This is in contrast to a computationally cheaper but less stable fully explicit scheme which uses forward differences at time $t_n$ (we will later consider the Crank-Nicolson method, which uses central differences at time $t_{n+\frac{1}{2}}$). Substituting the approximations

$$\frac{\partial V}{\partial S}(S_i, \tau_{n+1}) \approx \frac{V_{i+1}^{n+1} - V_{i-1}^{n+1}}{2\Delta S} \qquad (i=1,\dots,M-1, \quad n = 0, \dots, N-1),$$

$$\frac{\partial^2 V}{\partial S^2}(S_i, \tau_{n+1}) \approx \frac{V_{i+1}^{n+1} - 2 V_i^{n+1} + V_{i-1}^{n+1}}{(\Delta S)^2} \qquad (i=1,\dots,M-1, \quad n = 0, \dots, N-1),$$

$$ \frac{\partial S}{\partial \tau}(S_i, \tau_{n+1}) \approx \frac{V_i^{n+1} - V_i^n}{\Delta \tau} \qquad (i=0, \dots, M, \quad n = 0, \dots, N-1)$$

into our equation, we obtain for $i = 1,\dots, M-1$ and $n = 0, \dots, N-1$ the equation

$$ V_i^n = a_i V_{i-1}^{n+1} + b_i V_i^{n+1} + c_i V_{i+1}^{n+1} $$

where

$$a_i = -\frac{\Delta \tau}{2}(\sigma^2 i^2 - ri), \quad b_i = 1 + \Delta\tau(\sigma^2 i^2 + r), \quad c_i = -\frac{\Delta\tau}{2}(\sigma^2 i^2 + ri).$$


The initial condition $V(S,0) = (K-S)^+$ implies $V_i^0 = V(S_i, 0) = (K-S_i)^+$ for each $i = 0, \dots, M$, the boundary condition $V(0,\tau) = Ke^{-r\tau}$ implies $V_0^n = V(S_0, \tau_n) = V(0, \tau_n) = Ke^{-r\tau_n}$ for $n = 0, \dots N$ and the boundary condition $V(S_\text{max},\tau) = 0$ implies $V_M^n = V(S_\text{max},\tau_n) = 0$ for $n=0,\dots,N$. 


Thus for each $n=0,\dots,N-1$ we have the linear system 

$$
\begin{align}
\begin{pmatrix}
b_1 & c_1 & 0 & 0 & \cdots & 0 & 0 & 0 \\
a_2 & b_2 & c_2 & 0 & \cdots & 0 & 0 & 0\\
0 & a_3 & b_3 & c_3 & \cdots & 0& 0 & 0\\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
0 & 0 & 0 & 0 & \cdots & a_{M-2}& b_{M-2} & c_{M-2} \\
0 & 0 & 0 & 0 & \cdots & 0 & a_{M-1} & b_{M-1}
\end{pmatrix} 
\begin{pmatrix} V_1^{n+1} \\[4pt]
V_2^{n+1} \\
\vdots \\[4pt]
V_{M-2}^{n+1} \\[4pt]
V_{M-1}^{n+1}
\end{pmatrix} & = \begin{pmatrix} V_1^n - a_1 V_0^{n+1} \\[4pt] 
V_2^{n} \\
\vdots \\[4pt]
V_{M-2}^n \\[4pt]
V_{M-1}^n - c_{M-1}V_M^{n+1} 
\end{pmatrix} \nonumber \\[40pt]
& = \begin{pmatrix} V_1^n - a_1 Ke^{-r\tau_{n+1}} \\[4pt] 
V_2^{n} \\
\vdots \\[4pt]
V_{M-2}^n \\[4pt]
V_{M-1}^n
\end{pmatrix}\nonumber
\end{align}
$$

(which incorporates the two boundary conditions), along with the initial condition 

$$
\begin{pmatrix}
V_1^0 \\[4pt]
V_2^0 \\
\vdots \\[4pt]
V_{M-2}^0 \\[4pt]
V_{M-1}^0 
\end{pmatrix} = \begin{pmatrix}
(K-S_1)^+ \\[4pt]
(K-S_2)^+ \\
\vdots \\[4pt]
(K-S_{M-2})^+ \\[4pt]
(K-S_{M-1})^+ 
\end{pmatrix}. 
$$

Inverting the tridiagonal matrix above then allows us to iteratively find the value of $V_i^n$ at any time step $n$. 