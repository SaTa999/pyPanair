Disclaimer: This is a rough draft for `tutorial1.ipynb`. If you're willing to take the tutorial check the jupyter notebook instead.

# pyPanair Tutorial#1 Rectangular Wing  
---------------------------------------

In this tutorial we will perform an analysis of a rectangular wing with a NACA0012 airfoil.  
A brief overview of the procedure is listed below:  
1. Define the geometry of the wing using `wgs_creator.py`, and create input files `naca0012.wgs` and `naca0012.aux` for `panin`
2. Using the preprocessor `panin`, create an input file `a502.in` for `panair`
3. Run the analysis
4. Visualize the results from the analysis via `agps_converter.py`, `ffmf_converter.py`, and `calc_section_force.py`

## 1. Defining the geometry
---------------------------------------
First off, we will begin by defining the geometry of the rectangular wing.

The input geometry for `panair` (or more precisely, its preprocessor `panin`) is defined in the Langley Wireframe Geometry Standard (LaWGS) format. The format is described in [reference 1](https://ntrs.nasa.gov/search.jsp?R=19850014069).  

In a nutshell, the LaWGS format is a bundle of matrices, which are composed of m x n 3-dimensional points. Each matrix is referred to as a "network". Below is an example of a LaWGS file for a delta wing.

```sample1.wgs
puts 'The best way to log and share programmers knowledge.'
```

### References
---------------------------------------
[1] Craidon, C. B., "A Description of the Langley Wireframe Geometry Standard (LaWgs) Format," *NASA TM 85767*, 1985.
