Disclaimer: This is a rough draft for `tutorial3.ipynb`.
If you're willing to take the tutorial check the jupyter notebook instead.

# pyPanair Tutorial#3 Wing/body model
In this tutorial we will perform an analysis of a wing/body configuration.  

The model is based on the [AGARD-B model<sup>1</sup>](http://www.uwal.org/download/documents/agardcalmodelspecs.pdf)  
(The sting is not included, and the aft-body is slightly modified for this tutorial.)
![uppersurface]()

## 1.Defining the geometry

The model will be composed from eight networks:    

1. Wing  
2. Nose  
3. Mid-body
4. Aft-body
5. Body base
6. Wing wake
7. Body base wake
8. Body wake

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

Next, we'll define the tip chord. 
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
```python

```
```python

```



### Reference
1. Damljanović, D., Vitić, A., Vuković, Đ., and Isaković, J. 
"Testing of AGARD-B calibration model in the T-38 trisonic wind tunnel," *Scientific Technical Review*, Vol. 56, No. 2, 2006, pp. 52-62.