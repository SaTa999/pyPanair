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

```python

```

```python

```

```python

```

```python

```

```python

```



### Reference
1. Damljanović, D., Vitić, A., Vuković, Đ., and Isaković, J. 
"Testing of AGARD-B calibration model in the T-38 trisonic wind tunnel," *Scientific Technical Review*, Vol. 56, No. 2, 2006, pp. 52-62.