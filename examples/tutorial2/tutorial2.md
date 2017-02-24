Disclaimer: This is a rough draft for `tutorial2.ipynb`.
If you're willing to take the tutorial check the jupyter notebook instead.

# pyPanair Tutorial#2 Tapered Wing  
In this tutorial we will perform an analysis of a tapered wing.  
The wing is defined by five different wing sections at $\eta=0.000, 0.126, 0.400, 0.700, 1.000$.  

Below are the wing planform and airfoil stack, respectively.

![planform](planform.png)

```python
%matplotlib notebook
import matplotlib.pyplot as plt
from pyPanair.preprocess import wgs_creator
for eta in ("0000", "0126", "0400", "0700", "1000"):
    af = wgs_creator.read_airfoil("eta{}.csv".format(eta)) 
    plt.plot(af[:,0], af[:,2], "k-", lw=1.)
plt.plot((504.9,), (0,), "ro", label="Center of Gravity")
plt.legend()
plt.show()
```

(The wing is based on the [DLR-F4<sup>1</sup>](https://aiaa-dpw.larc.nasa.gov/Workshop1/files/agard-ar-303.pdf))

## 1.Defining the geometry
Just as we have done in tutorial 1, we will use the `wgs_creator` module to define the geometry of the wing.

First off, we initialize a `LaWGS` object.

```python
wgs = wgs_creator.LaWGS("tapered_wing")
```

Next, we create a `Line` object that defines the coordinates of the airfoil at the root of the wing.  
To do so, we will read a csv file that contains the coordinates of the airfoil, using the `read_airfoil` function.  

Three csv files, `root.csv`, `kink.csv`, and `tip.csv` have been prepared for this tutorial.  

Before we creating the `Line` object, we will take a quick view of these files.  
For example, the `root.csv` looks like ...

```python
import pandas as pd
pd.set_option("display.max_rows", 10)
pd.read_csv("root.csv")
```

The first and second columns `xup` and `zup` represent the x & z coordinates of the upper surface of the airfoil.  
The third and fourth columns `xlow` and `zlow` represent the x & z coordinates of the lower surface of the airfoil.  

The csv file must follow four rules:  
1. Data in the first row correspond to the x & z coordinates of the leading edge of the airfoil  
2. Data in the last row correspond to the x & z coordinates of the trailing edge of the airfoil  
3. For the first row, the coordinates `(xup, zup)` and `(xlow, zlow)` are the same  
4. For the last row, the coordinates `(xup, zup)` and `(xlow, zlow)` are the same (i.e. the airfoil has a sharp TE)  

### Reference
1. Redeker, G., "A selection of experimental test cases for the validation of CFD codes,"
 *AGARD AR-303*, 1994.
