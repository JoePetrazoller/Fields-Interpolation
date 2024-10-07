# Elastic fields interpolation

## Goal
The goal of this code is to interpolation the fields from irregular grid LAMMPS dump file, to a regular grid where the spatial step is way smaller.

## Interpolation

The interpolation is performed using the Gaussian Kernel :

$$\langle S\rangle=\dfrac{\Sigma_{i}\left(S_{i} \times W_{i}\right)}{\Sigma_{i} \: W_{i}}$$;

$$W_{i} = \dfrac{1}{\sqrt{2\pi\sigma ^2}} exp\left(\dfrac{-d_{i}^2}{2\sigma^2}\right)$$; 

$$d_{i} = \sqrt{\left( x_{1i} - x_{1} \right )^2 + \left( x_{2i} - x_{2} \right) ^2 }$$

Where $\langle S\rangle$ is the interpolated quantity, $W_{i}$ the weighting factors, $\sigma$ the standard deviation of the Gaussian law, ${x_1i}$ and ${x_2i}$ the coordinates of the atom i and $x_1$ and $x_2$ the coordinates of the considered node in the regular grid.

## How to use
1. Get a LAMMPS dump result
2. Place inside the same folder the dump file from LAMMPS. It should contain only the last increment.
