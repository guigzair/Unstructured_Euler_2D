
# Unstructured solver of the Euler equations
<!-- ABOUT THE PROJECT -->
## About The Project

The project is a a 2D solver of the Euler Equations on a unstructured mesh at the first order.
Given the Euler Equations:
$$\partial_t \begin{pmatrix}
\rho \\
\rho u \\
\rho E
\end{pmatrix} + \nabla .\begin{pmatrix}
\rho u \\
\rho u \otimes u + P I \\
 (\rho E + P)u
\end{pmatrix} = 0$$ 
With $E = e + \frac{1}{2}u^2$ and $e = \dfrac{P}{\rho (\gamma - 1)}$


Using the Rusanov Flux at each interface of each cells: 
$$F_{i,j} = \dfrac{F_i + F_j}{2} - max(C_i, C_j) (W_i - W_j)$$
Giving the updating using a Euler method for the integration : 
$$W_i^{t+dt} = W_i^{t} - \dfrac{dt}{A_i}\sum_{j\in N_j} F_{i,j}$$

The code is very basic and not optimized, boundaries are updated as outgoing supersonic for simplicity.
I made this code public because i couldn't find one on internet when I tried to do one by myself. 

![plot](./rho.jpg)

<!-- Requirements -->
## Requirements

numpy, meshpy, scipy, matplolib

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.


<!-- CONTACT -->
## Contact

Guillaume de Romemont - guillaume.romemont@protonmail.com





