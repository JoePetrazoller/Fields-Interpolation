# Elastic fields interpolation

## Goal
The goal of this code is to interpolation the Virial stress fields from irregular grid LAMMPS dump file, to a 2D regular grid where the spatial step is way smaller. Can be easily adapted to interpolate other quantities. The assumption that the stresses do not change along one direction direction is made (This direction is Z)

## Interpolation

The interpolation is performed using the Gaussian Kernel :

$$\langle S\rangle=\dfrac{\Sigma_{i}\left(S_{i} \times W_{i}\right)}{\Sigma_{i} \: W_{i}}$$

$$W_{i} = \dfrac{1}{\sqrt{2\pi\sigma ^2}} exp\left(\dfrac{-d_{i}^2}{2\sigma^2}\right)$$

$$d_{i} = \sqrt{\left( x_{1i} - x_{1} \right )^2 + \left( x_{2i} - x_{2} \right) ^2 }$$

Where $\langle S\rangle$ is the interpolated quantity, $W_{i}$ the weighting factors, $\sigma$ the standard deviation of the Gaussian law, ${x_1i}$ and ${x_2i}$ the coordinates of the atom i and $x_1$ and $x_2$ the coordinates of the considered node in the regular grid.

## How to use
1. Get a LAMMPS dump result, where the atom id, atom type (just for formatting, will not be used), x, y and z coordinates, the 6 component of the Virial stress tensor per atom and the voronoi volume per atom. The (0,0,0) coordinate should be at the bottom left of the atomistic box.
2. Place inside the same folder the dump file from LAMMPS. It should contain only the last increment.
3. Inside the interpolation file, set the sigma, the cutoff, the filename, and the size of the box in X and Y direction and launch the calculation.
4. For more efficiency, one can duplicate 6 times the code and run only one loop of the "for loop" on each file.
5. One .csv file is created per component. One can than easily open them using for example "contourf" function of Matlab.
