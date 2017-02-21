# pyPanair
A pre/post processor for Panair.  
"Panair" and "Panin" are not included in this repository.  
A copy of these software can be obtained at [Public Domain Aeronautical Software](http://www.pdas.com/contents15.html).  

## What can "pyPanair" do?  
List of things that are currently implemented in pyPanair (check the examples)
* Preprocessor that converts LaWGS format files to stl format (convert_wgs)  
* Postprocessor that converts agps format files to vtk (Legacy Paraview), vtm (Multiblock Paraview),
 and dat (multiblock Tecplot) format (convert agps)  
* Postprocessor that calculates the section force (local lift coefficient) from agps files (section_force_calculation)
* Postprocessor that converts ffmf files to csv format (convert_ffmf)  

List of things that will be implemented **soon<sup>TM</sup>**  
* Preprocessor for a simple wing geometry
* Preprocessor for a wing & body geometry

## Installation
Download the repository and type

```commandline
python setup.py install
```

or if you have git installed, simply type

```commandline
pip install git+https://github.com/SaTa999/pyPanair
```

## Requirements
pyPanair requires python 3  
(tests have only been performed for python 3.6)  
An [Anaconda3](https://www.continuum.io/) environment is recommended.

## Example
Example scripts and files are bundled in the "examples" directory.  
Run the scripts in each directory to check out how pyPanair works.   