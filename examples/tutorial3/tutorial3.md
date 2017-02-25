Disclaimer: This is a rough draft for `tutorial3.ipynb`.
If you're willing to take the tutorial check the jupyter notebook instead.

# pyPanair Tutorial#3 Wing/body model
In this tutorial we will perform an analysis of a wing/body configuration.  

The model is based on the [AGARD-B model<sup>1</sup>](http://www.uwal.org/download/documents/agardcalmodelspecs.pdf)  
(The sting is not included, and the aft-body is slightly modified for this tutorial.)
![uppersurface]()

## 1.Defining the geometry

The model will be composed from 12 networks:    

1. Wing  
2. Nose+Fore-body (upper surface)  
3. Nose+Fore-body (lower surface)  
4. Mid-body (upper surface)  
5. Mid-body (lower surface)  
6. Aft-body (upper surface)  
7. Aft-body (lower surface)  
8. Body base
9. Wing wake
10. Body base wake (upper)
11. Body base wake (lower)
12. Body wake

The geometry of the model is defined in terms of the body diameter `D`.  
For simplicity, it will be set at `D=1`.

### 1.1 Wing

First, initialize a `LaWGS` object.

```python
from pyPanair.preprocess import wgs_creator
wgs = wgs_creator.LaWGS("agardb_mod")
```

Next, we will create a `Line` object that defines the airfoil at the wing root.  
The airfoil of the AGARD-B model is a circular-arc airfoil of 4% thickness.  
A csv file has been prepared in advance. 
The coordinates in `circular_arc.csv` are normalized, so when using the `read_airfoil` function, we shall set the `expantion_ratio` at `2.598`.
By doing so, the x & z-axis coordinates of the resulting airfoil will be multiplied by `2.598`.

```python
root_airfoil = wgs_creator.read_airfoil("circular_arc.csv", span_pos=0.5, expansion_ratio=2.598)
%matplotlib notebook
import matplotlib.pyplot as plt
plt.plot(root_airfoil[:,0], root_airfoil[:,2], "s", mfc="None", mec="b")
plt.xlabel("x")
plt.ylabel("z")
plt.grid()
```

Since the apex of the nose is the origin, we'll need to shift the coordinates of the `root_airfoil`, so that the x-coordinate of its LE becomes `4.5`.
This will be done by with the `shift` method.

```python
```

Next, we define the tip chord. 
The tip of the wing is a point, but when defining it, we need to use the same number of points we used to define the wing root.
This means that we need to define a `Line` with `61` identical points.  

This can be easily done by using the `replace` method. By typing

```python
tip_airfoil = root_airfoil.replace(x=2.598, y=2.)
```

a copy of `root_airfoil` will be created. 
The xyz coordinates of each point in the copy will be replaced by `2.598`, `2.`, and `0.`, respectively.

The wing `Network` will be created and registered to `wgs` in the same way as tutorials 1 and 2.  

```python
wing = root_airfoil.linspace(tip_airfoil, num=30)
wgs.append_network("wing", wing, 1)

wing.plot_wireframe()
```

### 1.2 Nose+Fore-body

The nose profile is described by
$$r=\frac{x}{3}[1-\frac{1}{9}(\frac{x}{D})^2+\frac{1}{54}(\frac{x}{D})^3]$$
$r$ is the radius and $x$ is the x-coordinate.

The fore-body is a cylinder with a radius of $D$.

First, we define the `Line` for the nose and fore-body at `z=0.`.

```python
import numpy as np

n_nose = 15
x_nose = np.linspace(0., 3., num=n_nose)
y_nose = x_nose/3 * (1 - 1/9*(x_nose)**2 + 1/54*(x_nose)**3)
nose_line = wgs_creator.Line(np.zeros((n_nose, 3)))
nose_line[:,0] = x_nose
nose_line[:,1] = y_nose

fbody_p1 = wgs_creator.Point(nose_line[-1])
fbody_p2 = fbody_p1.replace(x=4.5)
fbody_line = fbody_p1.linspace(fbody_p2, num=5)

nose_fbody_line = nose_line.concat(fbody_line)
```
Other `Lines` that define the nose and fore-body will be created by rotating `fbody_line`.
The upper/lower surfaces can be created by typing

```python
nose_fbody_up = list()
for i in np.linspace(0, 90, 7):
    line = nose_fbody_line.rotx(rotcenter=nose_fbody_line[0], angle=i)
    nose_fbody_up.append(line)
nose_fbody_up = wgs_creator.Network(nose_fbody_up)

nose_fbody_low = list()
for i in np.linspace(-90, 0, 7):
    line = nose_fbody_line.rotx(rotcenter=nose_fbody_line[0], angle=i)
    nose_fbody_low.append(line)
nose_fbody_low = wgs_creator.Network(nose_fbody_low)

wgs.append_network("n_fb_up", nose_fbody_up, 1)
wgs.append_network("n_fb_low", nose_fbody_low, 1)
```

### 1.3 Mid-body

The mid-body is a cylinder with a hole at the section where the body and wing intersect.  

First, we define the upper surface.

```python
wingroot_up, wingroot_low = root_airfoil.split_half()
mbody_line = wingroot_up.replace(z=0.)

mbody_up = list()
for i in np.linspace(90, 0, 7)[:-1]:
    line = mbody_line.rotx((0,0,0), angle=i)
    mbody_up.append(line)
mbody_up.append(wingroot_up)
mbody_up = wgs_creator.Network(mbody_up)
wgs.append_network("mbody_up", mbody_up, 1)
```

Next, we define the lower surface.

```python
wingroot_low = wingroot_low.flip()
mbody_low = list()
mbody_low.append(wingroot_low)
for i in np.linspace(0, -90, 7)[1:]:
    line = mbody_line.rotx((0,0,0), angle=i)
    mbody_low.append(line)
mbody_low = wgs_creator.Network(mbody_low)
wgs.append_network("mbody_low", mbody_low, 1)
```

### 1.4 Aft-body
  
In the original AGARD-B model, the aft-body is a cylinder.  
However, in this tutorial it will be a circular truncated cone, for the sake of learning how to attach body wakes.  

First, we define the `Line` for the aft-body at `z=0.`.

```python
aft_body_p1 = wgs_creator.Point(root_airfoil[0])
aft_body_p2 = aft_body_p1.shift((1.402, -0.05, 0.))
aft_body_line = aft_body_p1.linspace(aft_body_p2, num=5)
```

After that, we define the `Networks` for the upper and lower surfaces.

```python
aft_body_up = list()
for i in np.linspace(0, 90, num=7):
    line = aft_body_line.rotx((0,0,0), angle=i)
    aft_body_up.append(line)
aft_body_up = wgs_creator.Network(aft_body_up)
wgs.append_network("abody_up", aft_body_up, 1)

aft_body_low = list()
for i in np.linspace(-90, 0, num=7):
    line = aft_body_line.rotx((0,0,0), angle=i)
    aft_body_low.append(line)
aft_body_low = wgs_creator.Network(aft_body_low)
wgs.append_network("abody_low", aft_body_low, 1)
```

### 1.5 Body base

The body base is a circle.
The arc is defined by `edge3` from the `Networks` `aft_body_up` and `aft_body_low`  
We'll use the `edge` method to create `Line` objects from these `Networks`.

```python
body_base_line_up = aft_body_up.edge(3)
body_base_line_up = body_base_line_up.flip()
body_base_line_low = aft_body_low.edge(3)
body_base_line_low = body_base_line_low.flip()
body_base_line = body_base_line_up.concat(body_base_line_low)
body_base_line2 = body_base_line.replace(x=8.5, y=0., z=0.)
body_base = body_base_line.linspace(body_base_line2, num=3)
```

Note that when registering a body base, the boundary type should be `5`.

### 1.6 Wakes

For this model we need to attach 3 types of wakes: 
1. Wing wake
2. Body base wake
3. Body wake

The wing wake can be created in the same way as tutorials 1 and 2.  

```python
wake_length = 2.3093 * 50
wing_wake = wing.make_wake(3, wake_length)
```

The body base wake will be attached to `edge3` of the `Networks` `aft_body_up` and `aft_body_low`.

```python
body_base_wake_up = aft_body_up.make_wake(3, wake_length)
wgs.append_network("b_wake_up", body_base_wake_up, 18)
body_base_wake_low = aft_body_up.make_wake(3, wake_length)
wgs.append_network("b_wake_low", body_base_wake_low, 18)
```

The body wake is a wake that connects the wing wake and body base wake. 
We will attach it to `edge4` of `aft_body_up`.  

```python
body_wake = aft_body_up.make_wake(4, wake_length)
```

Note that corner `1` of the body base wake must lie at the intersection of the wing TE and the inboard edge of the wing wake.
Also, the boundary type for a body wake is `20`.

(Read 3.3.2.5 of [Vol.2 of the user manual<sup>2</sup>](https://ntrs.nasa.gov/search.jsp?R=19920013622) for more information on body wakes.)  

Next, we create a stl file and check for errors.  

```python
wgs.create_stl("agardb_mod.stl", include_wake=True)
```

Finally, we create input files for `panin`.

```python
wgs.create_wgs("agardb_mod.wgs")
wgs.create_aux("agardb_mod.aux", alpha=5., mach=1.4, cbar=2.3093, span=4., sref=6.928, xref=5.943, zref=0.)
```

## 2. Analysis

After creating the geometry of the model, run the analysis.  
(If you are unaware of how to do so, check tutorial 1 or 2.)

## 3. Validation and visualization

First, we'll validate the results by comparing the aerodynamic coefficients to reference 1.  

```python
from pyPanair.postprocess import read_ffmf
read_ffmf()
```

According to reference 1, the aerodynamic coefficients for an AGARD-B model at an angle of attack of 5.0 degrees, mach number of 5.0, and a Reynolds number of around 10.4 million is,  
$$\begin{align}C_L&=0.224 \\
C_D&=0.0435 \\
C_M&=-0.019\end{align}$$  
The lift and drag coefficients are fairly close to the results from reference 1, considering the fact that we've slightly tampered with the shape of aft-body.  
On the other hand, the pitching moment coefficient does not match with reference 1. 
This may be due to either the modified aft-body or the coarse paneling.

Visualization can be done in the same way as tutorials 1 and 2.  

```python
from pyPanair.postprocess import write_vtk, write_tec
write_vtk(n_wake=4)
write_tec()

import pandas as pd
from pyPanair.postprocess import calc_section_force
calc_section_force(aoa=5., mac=2.3093, rot_center=(5.943,0,0), casenum=1, networknum=1)

section_force = pd.read_csv("section_force.csv")
section_force

section_force = section_force.dropna()
plt.plot(section_force.pos / 2, section_force.cl * section_force.chord, "s", mfc="None", mec="b")
plt.xlabel("spanwise position [normalized]")
plt.ylabel("$C_l$ * chord")
plt.grid()
plt.show()
```

This is the end of tutorial 3.  

### Reference
1. Damljanović, D., Vitić, A., Vuković, Đ., and Isaković, J. 
"Testing of AGARD-B calibration model in the T-38 trisonic wind tunnel," *Scientific Technical Review*, Vol. 56, No. 2, 2006, pp. 52-62.
2. Sidwell, K. W., Baruah, P. K., Bussoletti, J. E., Medan, R. T., Conner, R. S., and Purdon, D. J.,
"PAN AIR: A computer program for predicting subsonic or supersonic linear potential flows about arbitrary configurations using a higher order panel method. Volume 2: User's manual (version 3.0),"
 *NASA CR 3252*, 1990.