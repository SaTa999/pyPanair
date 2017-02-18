# pyPanair
A pre / post processor for PANAIR.  
"a502" (PANAIR), "panin", and "makewgs" are not included in this repo.  
A copy of each of these software can be obtained at [Public Domain Aeronautical Software](http://www.pdas.com/contents15.html).  

## What can "pyPanair" do?  
Currently, post processors for agps files (agps_converter) and ffmf files (ffmf_converter) are implemented.  

List of things that will be implemented **soon<sup>TM</sup>**  
* Preprocessor for a simple wing geometry
* Preprocessor for wing & body geometry
* Postprocessor for calculating the spanwise local lift distribution

## Installation
Download the repository and type

```commandline
python setup.py install
```

or if you have git installed, type

```commandline
pip install git+https://github.com/SaTa999/pyPanair
```

## Requirements
pyPanair requires python 3  
(tests have only been performed for python 3.6)  
An [Anaconda3](https://www.continuum.io/) environment is recommended.
