Disclaimer: This is a rough draft for `tutorial1.ipynb`. If you're willing to take the tutorial check the jupyter notebook instead.

# pyPanair Tutorial#1 Rectangular Wing  
In this tutorial we will perform an analysis of a rectangular wing with a NACA0012 airfoil.
A brief overview of the procedure is listed below:  

1. Define the geometry of the wing using `wgs_creator.py`, and create input files `naca0012.wgs` and `naca0012.aux` for `panin`
2. Using the preprocessor `panin`, create an input file `a502.in` for `panair`
3. Run the analysis
4. Visualize the results from the analysis via `agps_converter.py`, `ffmf_converter.py`, and `calc_section_force.py`

## 1. Defining the geometry
### 1.1 LaWGS Format
First off, we will begin by defining the geometry of the rectangular wing.

The input geometry for `panair` (or more precisely, its preprocessor `panin`) is defined in the Langley Wireframe Geometry Standard (LaWGS) format. The format is described in [reference 1](https://ntrs.nasa.gov/search.jsp?R=19850014069).  

In a nutshell, LaWGS files are a bundle of 3 dimensional arrays, which are referred to as "networks". A network is composed of "lines", which in turn is composed of 3-dimensional "points". If a network has `m` lines, and each line of has `n` points, the shape of the corresponding 3d array will be `(m, n, 3)`. Below is an example of a LaWGS file for a delta wing.

```sample1.wgs
deltawing created from wgs_creator
wing
1 3 5 0   0 0 0   0 0 0    1 1 1  0
1.0000000e+01 0.0000000e+00 0.0000000e+00
5.0000000e+00 0.0000000e+00 1.0000000e+00
0.0000000e+00 0.0000000e+00 0.0000000e+00
5.0000000e+00 0.0000000e+00 -1.0000000e+00
1.0000000e+01 0.0000000e+00 0.0000000e+00
7.5000000e+00 1.0000000e+01 0.0000000e+00
5.0000000e+00 1.0000000e+01 5.0000000e-01
2.5000000e+00 1.0000000e+01 0.0000000e+00
5.0000000e+00 1.0000000e+01 -5.0000000e-01
7.5000000e+00 1.0000000e+01 0.0000000e+00
5.0000000e+00 2.0000000e+01 0.0000000e+00
5.0000000e+00 2.0000000e+01 0.0000000e+00
5.0000000e+00 2.0000000e+01 0.0000000e+00
5.0000000e+00 2.0000000e+01 0.0000000e+00
5.0000000e+00 2.0000000e+01 0.0000000e+00
```

The first row displays the title of the LaWGS file.
```
deltawing created from wgs_creator
```

The second row displays the name of the network.
```
wing
```

The third row lists the parameters of the network.
```
1 3 5 0   0 0 0   0 0 0    1 1 1  0
```
The definitions of the first three numbers are as follows:
* "1": the id of the network
* "3": the number of lines in the network
* "5": the number of points in each line

The remaining 11 numbers, `0   0 0 0   0 0 0    1 1 1  0`, define the local and global axes. (Read reference 1 for more information.)  

The fourth and subsequent lines, define the coordinates of each point. For example, the fourth line, `1.0000000e+01 0.0000000e+00 0.0000000e+00`, means that the coordinate of the first point is `x=1., y=0., z=0.`.

The wireframe defined by the above file looks like this ...

```python
%matplotlib notebook
from pyPanair.preprocess import wgs_creator as wg
delta_wing = wg.read_wgs("sample1.wgs")
delta_wing._networks[0].plot_wireframe(show_normvec=False, show_corners=False, show_edges=False)
```

### 1.2 Creating a LaWGS file using pyPanair
Next, we will create a LaWGS file using the `wgs_creator` module of `pyPanair`. In the `wgs_creator` module, LaWGS files are created using objects that derive from the four classes, `LaWGS`, `Network`, `Line`, and `Point`. A brief explanation of these classes are written below.  
* `LaWGS`: A class that represents a LaWGS format geometry as a list of `Networks`. Can be used to read/write LaWGS files.
* `Network`: A class that describes a network as a 3-dimensional array.
* `Line`: A class that represents a line as a 2-dimensional array.
* `Point`: A class that defines the xyz coordinates of a point. A 1-dimensional array is used to define the coordinates.

Now we will begin the actual work of creating a LaWGS file. First, we start of by initializing a `LaWGS` object. The title of the LaWGS object will be `"NACA0012"`.

```python
wgs = wgs_creator.LaWGS("NACA0012")
```

In the next step, the `Network` of the wing will be defined by interpolating two `Lines`, the `root_airfoil` and `tip_airfoil`. The `root_airfoil`, which is an NACA0012 airfoil, can be constructed using the `naca4digit` method.

```python
root_airfoil = wgs_creator.naca4digit("0012", num=25, chord=100., span_pos=0.)
```

The resulting airfoil has a chord length of `100.`, and its spanwise position (i.e. y-axis coordinate) is `0.`.
The upper and lower surfaces are each represented with `25` points.
The x and z coordinates of the airfoil look like below:

```python
import matplotlib.pyplot as plt
plt.plot(root_airfoil[:,0], root_airfoil[:,2], "s", mfc="None", mec="b")
plt.xlabel("x")
plt.ylabel("z")
plt.grid()
```

The `tip_airfoil` can be defined in the same manner as the `root_airfoil`.
However, this time, the spanwise position will be `span_pos=300.`.

```python
tip_airfoil = wgs_creator.naca4digit("0012", num=25, chord=100., span_pos=300.)
```

The `Network` of the wing will be defined by interpolating to the two `Lines`, `root_airfoil` and `tip_airfoil`.
To do so, we simply use the `linspace` method.

```python
wing = root_airfoil.linspace(tip_airfoil, num=20)
```

The wing `Network` will have `20` lines, which are linear interpolations of the `root_airfoil` and `tip_airfoil`.
`Networks` can be visualized using the `plot_wireframe` method.

```python
wing.plot_wireframe()
```

Along with coordinates of each point in the `Network`, the corners (e.g. `1`) and edges (e.g. `edge1`) are displayed.
(Read reference 1 for more information on network corners and edges.)
Also, an arrow indicating the front side of the network is depicted.
(Details of "front side" will be mentioned later.)

After defining the `Network` for the wing, we register it to `wgs` using the `append_network` method.

```python
wgs.append_network("wing", wing, 1)
```

The first variable, `"wing"`, is the name of the network.
The second variable, `wing`, is the `Network` we are registering.
The third variable, `1`, is the boundary type of the network. If the network represents a solid wall, the type is `1`.
(Read [reference 2](https://docs.google.com/file/d/0B2UKsBO-ZMVgS1k5VElNamx1cUk/edit) for more information on the different types of boundary conditions.)

The next process will be to define the geometry of the wingtip.
To do so, we will split the `tip_airfoil` into upper and lower halves, and linearly interpolate them.
All of this can be done by typing

```python
wingtip_upper, wingtip_lower = tip_airfoil.split_half()
wingtip_lower = wingtip_lower.flip()
wingtip = wingtip_upper.linspace(wingtip_lower, num = 5)
```

The wing tip will look like ...

```python
wingtip.plot_wireframe()
```

The `wingtip` will also be registered to `wgs`.

```python
wgs.append_network("wingtip", wingtip, 1)
```

In addition to the `wing` and `wingtip`, we also must define the "wake" of the wing.
In `Panair`, the wake is defined as a square network stretching out from the trailing edge of the wing.
The length of the wake should be about 25 to 50 times the length of the reference chord of the wing.
The wake can be easily defined by using the method `make_wake`.

```python
wingwake = wing.make_wake(edge_number=3, wake_length=50*100.)
```

The `edge_number` variable, means that we are attaching the wake to `edge3` of the `Network` `wing`.
The `wake_length` variable defines how long the wake is.
In this case, the wake stretches from `x=100.` (the TE of the wing) to `x=5000.`.

The `wingwake` will also be registered to `wgs`.

```python
wgs.append_network("wingwake", wingwake, 18)
```

Notice that this time we are setting the boundary type as `18`.
Boundary type `18` is used to define the that network is a "wake" emitted from sharp edges.

Now that we have finished defining the geometry of the rectangular wing,
we will check to see if there are any errors in the model we've constructed.

To make this task easier, we will write the geometry into a STL (STereoLithography) format.
To do so, type

```python
wgs.create_stl("naca0012.stl")
```

A stl file named `naca0012.stl` should be created in the current working directory.
Open this file with a stl viewer. (I recommend [Materialise MiniMagics 3.0](http://www.materialise.co.jp/minimagics-stl-viewer-0))

Below is a screen shot of the stl.

![stl_screenshot](stl_screenshot.png)

Using the stl viewer, we should watch out for the following four points:

1. The stl model is watertight (i.e. No holes between each network)  
(Note that there seems to be a hole at the root of the wing.
 However, this model is symmetric with respect to the x-z plane, so the hole actually doesn't exist.)
2. There are no intersecting networks (i.e. No networks crossing over each other)
3. The front side of the network (excluding the wake network) is facing outwards.  
In the picture above, the grey side of the stl (which correspond to the front side of the networks) is facing the flow.
4. The corners of adjacent networks match each other

If you've followed the instructions, the `naca0012.stl` should fulfill all four points. (If not try again.)

Finally, we will write the input files for `panin`.
This can be accomplished, using the methods, `write_wgs` and `write_aux`.

```python
wgs.create_wgs("naca0012.wgs")
wgs.create_aux("naca0012.aux", alpha=2, mach=0.2, cbar=100., span=600.,
               sref=60000., xref=25., zref=0.)
```

Two files, `naca0012.wgs` and `naca0012.aux` should be crated in the current directory.  
`naca0012.wgs` is a file that defines the geometry of the model.  
`naca0012.aux` is a file that defines the analysis conditions.  
The definition of each variable is listed below: 
* `alpha`: The angle of attack (AoA) of the analysis.  
When conducting multiple cases, the variable should be a tuple defining the AoAs. (e.g. `(2, 4, 6, 8)`) 
Up to four cases can be conducted in one run.
* `mach`: The mach number of the flow
* `cbar`: The reference chord of the wing
* `span`: The reference span of the wing
* `sref`: The reference area of the wing
* `xref`: The x-axis coordinate of the center of rotation
* `zref`: The z-axis coordinate of the center of rotation

## 2. Creating an input file for `panair`

In this chapter we will create an input file for `panair`, using its preprocessor `panin`.  

If you do not have a copy of `panin`, download it from [PDAS](http://www.pdas.com/panin.html), and compile it.  
(I'm using cygwin, so the code blocks below are for cygwin environment.)

```bash
gfortran -o panin.exe panin.f90
```

After compiling `panin`, place `panin.exe`, `naca0012.wgs`, and `naca0012.aux` under the `tutorial1/panin/` directory.  
Then run `panin`.

```bash
./panin
```

It will display 
```bash
 Prepare input for PanAir
  Version 1.0 (4Jan2000)
 Ralph L. Carmichael, Public Domain Aeronautical Software
 Enter the name of the auxiliary file: 
```

Enter `naca0012.aux`. If everything goes fine, it should display
```bash
          10  records copied from auxiliary file.
           9  records in the internal data file.
  Geometry data to be read from NACA0012.wgs                                                                    
 Reading WGS file...
 Reading network wing
 Reading network wingtip
 Reading network wingwake
 Reading input file instructions...
 Command  1 MACH 0.2
 Command 11 ALPHA 2
 Command  6 cbar 100.0
 Command  7 span 600.0
 Command  2 sref 60000.0
 Command  3 xref 25.0
 Command  5 zref 0.0
 Command 35 BOUN 1 1 18
 Writing PanAir input file...
  Files a502.in added to your directory.
 Also, file panin.dbg
 Normal termination of panin, version 1.0 (4Jan2000)
 Normal termination of panin
```

and `a502.in` (an input file for `panair`) should be created under the current directory.


## 3. Running `panair`

Now its time to run the analysis.  
If you do not have a copy of `panair`, download it from [PDAS](http://www.pdas.com/panair.html), and compile it.
(A bunch of warnings will appear, but it should work.)

```bash
gfortran -o panair.exe -Ofast -march=native panair.f90
```

Place `panair.exe` and `a502.in` under the `tutorial1/panair/` directory, and run `panair`.

```bash
./panair
```

`panair` will display the below text.

```bash
 Panair High Order Panel Code, Version 15.0 (10 December 2009)
 Enter name of input file:
```

Enter `a502.in`. The analysis will end in a few seconds.  
After the analysis ends, the output files such as , `panair.out`, `agps`, and `ffmf` will be created in the current directory.  

`panair.out` contains the output of the whole analysis (e.g. source and doublet strength of each panel)  
`agps` contains the surface pressure distribution for each case  
`ffmf` contains the aerodynamic coefficients for each case  

* Warning: Along with the output files, you will also notice the existence of intermediate files (e.g. `rwms01`).  
Users should always delete these intermediate files when running new cases.  
(To do so, run `clean502.bat` or `clean502.sh` which is contained in the archive file `panair.zip`)

## 4. Validating and visualizing the output

### 4.1 Validation
In this chapter we will visualize the results of the analysis, but before we do so, we will validate the results by checking the aerodynamic coefficients.  

Open the `ffmf` file contained in the `tutorial1/panair/` directory with a text editor.  
After the headers of the file, you shall see 

```
 sol-no     alpha      beta            cl           cdi            cy            fx            fy            fz
                                                                                 mx            my            mz            area
 ------   -------   -------       -------       -------     ---------     ---------     ---------     ---------    ------------

      1    2.0000    0.0000       0.16011       0.00068       0.00000      -0.00491       0.00000       0.16004
                                                                            0.00000       0.00067       0.00000    123954.39335
```

This area shows the aerodynamic coefficients of each case.  
A brief explanation of each column is listed below:  

* `sol-no`: The case number
* `alpha`: The AoA of the case
* `beta`: The side slip angle of the case
* `cl`: The lift coefficient of the entire geometry
* `cdi`: The induced drag coefficient of the entire geometry
* `cy`: The side force coefficient of the entire geometry
* `fx, fy, fz`: The non-dimensional force in x, y, z direction, respectively
* `mx, my, mz`: The non-dimensional torque in x, y, z direction, respectively

Notice that `cl` and `fz` do not match.
This is because, `cl` is obtained from a trefz plane analysis, whereas `fz` is obtained by integrating the surface pressure of the geometry.
The same can be said for `cdi` and `fx`.

According to the lifiting line theory<sup>(3</sup>, when the AoA is $\alpha\mathrm{[rad]}$, the lift coefficient ($C_L$) and induced drag coefficient ($C_{D_i}$) for a untwisted uncambered rectangular wing with an aspect ratio of $6$ is 
$$C_L = 0.9160\frac{\pi^2}{2}\alpha$$
$$C_{D_i} = 0.8744\frac{\pi^3}{24}\alpha^2$$
In our analysis, the AoA is $0.1047 \mathrm{[rad]}$, so the lift and drag coefficients should be $C_L = 0.4734$ and $C_{D_i} = 0.01239$.  
The analysis predicted a fairly close value of $C_L = 0.4785$ and $C_{D_i} = 0.01266$.  

### 4.2 Visualization of the surface pressure distribution
Now we shall move on to the visualization of the result.  

First, we begin by converting the `agps` file into a format that can be used in common visualization applications.  
The `agps` file can be converted into three formats:
  
1. `vtk`: Legacy paraview format  
2. `vtu`: Multi-block paraview format  
3. `dat`: Multi-block tecplot format

In this tutorial we choose the `vtk` format.
To convert the `agps` file, first move the `agps` file to the `tutorial1/` directory.  
Then, use the `write_vtk` method of `pyPanair`. (If you wish to use tecplot, enter `write_tec` instead of `write_vtk`.)  

```python
from pyPanair.postprocess import agps_converter
agps_converter.write_vtk(n_wake=1)
```

(The `n_wake` variable is used to input the number of wakes. 
For example, if we enter `1`, the last `1` networks in the geometry will not be included in the output `vtk` file.)

`agps.vtk` should be created in the current directory.  
It can be open using [ParaView](http://www.paraview.org/).  

Below is a screen shot of ParaView.  

### 4.3


### References
---------------------------------------
1. Craidon, C. B., "A Description of the Langley Wireframe Geometry Standard (LaWgs) Format," *NASA TM 85767*, 1985.
2. Saaris, G. R., "A502I User's Guide-PAN AIR Technology Program for Solving Potential Flow about Arbitrary Configurations," 1992.
3. Moran, J., *An Introduction to Theoretical and Computational Aerodynamics*, John Wiley & Sons, Inc., 1984.

