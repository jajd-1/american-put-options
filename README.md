# American put option pricing as a parabolic free boundary problem 

The goal of this project is to formulate the American put option pricing problem mathematically as a parabolic obstacle/free boundary problem, solve the problem numerically and explore how the free boundary is affected by parameters. 

More precisely, we first construct and validate our numerical solution produced by the finite difference method (FDM) against the closed-form solution to the Black-Scholes PDE in the case of European put options. We then incorporate the obstacle imposed in the American case, and re-solve the problem numerically to extract the free boundary, which determines when it is optimal for the holder to exercise the option. In the process we address numerical challenges occuring near the free boundary. We then study how the free boundary changes as parameters such as the volatility, interest rate and maturity date are varied. 

## European put options

A European put option gives the holder the right (but not the obligation) to sell an asset at a specified strike price $K$ at a specified maturity time $T$. If $S_t$ denotes the price of the asset at time $t \leq T$, then the payoff at maturity is 

$$(K-S_T)^+.$$

The price of the option at times $t < T$ may be considered as a function of $S$ and $t$ and, under certain assumptions, is given by the solution to the Black-Scholes equation 

$$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS\frac{\partial V}{\partial S} - rV = 0$$

with terminal condition $V(S,T) = (K-S)^+$. Here, $\sigma$ denotes the volatility and $r$ is the risk-free rate. We also have the boundary conditions $\lim_{S\rightarrow \infty} V(S,t) = 0$ (the put option is essentially worthless when the underlying price is far above the strike price) and $V(0,t) = K e^{-r(T-t)}$ (if the stock is worthless then the put option is worth its strike price discounted from maturity back to time $t$). 

Making the substitutions

$$\tau = T - t, \quad u = Se^{r\tau} \quad \text{and}\quad x = \ln\bigg(\frac{S}{K}\bigg) + \bigg(r - \frac{1}{2}\sigma^2\bigg)\tau$$

transforms the Black-Scholes equation into the diffusion equation

$$\frac{\partial u}{\partial \tau} = \frac{1}{2}\sigma^2 \frac{\partial^2 u}{\partial x^2},$$

which can be solved explicitly and transformed back to the original variables. 

**To do:** add explicit formula.

## American put options

An American put option gives the holder the additional flexibility to exercise the option at any time before maturity. The first change required in the PDE formulation is the boundary condition at $S=0$: since one no longer has to wait until maturity, an American put option is simply worth its strike price if the stock price hits zero, i.e. $V(0,t) = K$. The second change is more interesting: rather than just the terminal condition $V(S,T) = (K-S)^+$, the option price must also satisfy

$$V(S,t) \geq (K-S)^+.$$

at all times $t \leq T$. Indeed, if this weren't the case then there would be an arbitrage opportunity: one could buy the asset at price $S$, purchase the put option at price $V$ and immediately exercise it to sell the asset at price $K$, resulting in a profit of $K - S - V > 0$.

Therefore, the price of an American put option is characterised by the constraint $V(S,t) \geq (K-S)^+$ together with the condition that the Black-Scholes equation is satisfied in the region where $V(S,t) > (K-S)^+$, which we refer to as the *continuation region*. In the *exercise region* where $V(S,t) = (K-S)^+$, it is optimal for the holder to exercise the option. The interface between these two regions is the so-called free boundary, which determines for each time $t$ the threshold asset price below which exercise is optimal and above which the Black-Scholes equation governs the option price. 

We can picture the free boundary as a curve $\{(S,t): S\in[0,\infty), t\in[0,T]\}\subset\mathbb{R}^2$, although stricty speaking this requires one to prove something about the structure and regularity of free boundaries for Lipschitz obstacles. We put aside the analytic aspects of this free boundary problem for now, and consider instead numerical approaches to the problem. 

## Numerical analysis 


