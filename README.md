# model H







For simplicity, we first look at the continuum equations of motion for
($\phi$, $\pi$) without dissipation or noise:
$$\begin{aligned}
    \dot{\phi} &= - \nabla\_\mu \left( \phi \pi\_{T,\mu} \right) = - \pi\_{T,\mu}  \nabla\_\mu \phi, \\
    %
    \dot{\pi}\_{T,i} &= - \hat{{\cal T}}\_{ij} \nabla_k \left( \pi\_{T,k} \pi\_{T,j} + \nabla_k \phi \nabla_j \phi \right) \\
    %
    & = - \hat{{\cal T}}\_{ij} \left\[ \nabla_k \left(\pi\_{T,k} \pi\_{T,j} \right) + \nabla_j \phi \nabla^2 \phi \right\] 
\end{aligned}$$

The time derivative of the Hamiltonian is (assumed *ρ* = 1),
$$\begin{aligned}
    \dot{{\cal H}} = \int d^3x \\ \left\[ - \dot{\phi} \\ \nabla^2 \phi + \pi\_{T,i} \dot{\pi}\_{T,i} + V'(\phi) \\ \dot{\phi} \right\]
\end{aligned}$$


$$\begin{aligned}
    \phi^{n+1/3}\_\mu &= \phi^n\_\mu + \Delta t \\ {\cal F}\_\mu (\phi^n\_\mu), \\
    %
    \phi^{n+2/3}\_\mu &= \frac{3}{4} \\ \phi^n\_\mu + \frac{1}{4} \\ \phi^{n+1/3}\_\mu + \frac{\Delta t}{4} \\ {\cal F}\_\mu\left( \phi^{n+1/3}\_\mu \right), \\
    %
    \phi^{n+1}\_\mu &= \frac{1}{3} \\ \phi^{n}\_\mu + \frac{2}{3} \\ \phi^{n+2/3}\_\mu + \frac{2}{3} \\ \Delta t \\ {\cal F}\_\mu \left( \phi^{n+2/3}\_\mu \right).
\end{aligned}$$
