## **Linking Number in Python**

Function to compute the Linking number between two 3-dimentional oriented curves in python.

<img src="https://github.com/gon-uri/linking_number/blob/master/fig.PNG" alt="figure" width="550">

---
#### **Brief description**
The function first interpolate both curves using the a B-spline method from `scipy` library. Then it projects both curves to a 2-dimentional plane and computes the instersection of the curves in that plane. For each intersection, it computes the vector product of the projected oriented curves on that point to check wether it contributes `+` or `-`. Finally it compute the sum all the contributions to obatain the Linking Number.
For a detailed description of what a Linking Number is and how it is computed, please refer to out paper [here](https://aip.scitation.org/doi/10.1063/5.0013714).

---
#### **Main Function**

The parameters to the main function `linking_number` are:

```
    Inputs
    curve_1: array_like
        array of shape (3, N1) with N1 the number of points in curve 1
    curve_2: array_like
        array of shape (3, N2) with N2 the number of points in curve 2
    projection: {'XY','ZX','YZ','AUTO'}, optional
        Projection plane where the intersections will be computed (default 'AUTO')
    puntos_curva: int, optional
        Number of interpolation for B-spline method (default 5000)
    margin:  int, optional
        Distance from intercetption to compute the vector (default 10)
    verbose: boolean, optional
        Print information about the computation (default False)

    Outputs
    total: float
        Linking Number (it should always be an integer, if not, check changing parameters)
    coords_1: 
        Coordinates of intesection points in the fist dimention of the projection plane (XY','ZX' or'YZ').
    coords_2: 
        Coordinates of intesection points in the second dimention of the projection plane (XY','ZX' or'YZ').
```
---
#### **Citing**
If you use this function in your work, please remember to cite it using:

```
G. Uribarri, and G. B. Mindlin, "The structure of reconstructed flows in latent spaces", Chaos 30, 093109 (2020)

```
Bibtex:
```
@article{doi:10.1063/5.0013714,
author = {Uribarri,Gonzalo  and Mindlin,Gabriel B. },
title = {The structure of reconstructed flows in latent spaces},
journal = {Chaos: An Interdisciplinary Journal of Nonlinear Science},
volume = {30},
number = {9},
pages = {093109},
year = {2020},
doi = {10.1063/5.0013714},
URL = { 
        https://doi.org/10.1063/5.0013714
},
eprint = { 
        https://doi.org/10.1063/5.0013714  
}
}
```

