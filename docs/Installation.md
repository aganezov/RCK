# Installation


### Contents: 
* [Virtual environment](#virtual-environment)
* [conda](#conda-recommended)
* [pip](#pip)
* [source](#source)
* [Gurobi](#gurobi)
    * [Python bindings for Gurobi](#python-bindings-for-gurobi)
* [Executables](#executables)

### virtual environment
We recommend that RCK is installed in the isolated virtual environments.
Virtual environments can be created via `anaconda` and `python` (w.r.t. RCK, which is written in Python).

To create a virtual environment (named `rck-env`) with anaconda, run the following command:
````bash
conda create --name rck-env python=3.7
```` 

To create a virtual environment (named `rck-env`) with python, run the following command:
````bash
python -m venv rck-env
````

If virtual environments are used (which, again, we recommend), we assume that the environment is activated.

### conda (recommended)
Run the following conda command, to install RCK:
````bash
conda install -c aganezov rck
````

Installation via conda automatically takes care of Gurobi python bindings (refer to respective [subsection](#python-bindings-for-gurobi)), and everything shall work from this part (assuming that Gurobi is correctly installed and working).

### pip

Run the following command, to install RCK:
````bash
pip install rck
````

**WARNING**: this installation does take care of python bindings for Gurobi. Please, refer to respective [subsection](#python-bindings-for-gurobi) on how that can be addressed.

### source

First, download the source code. Example is shown below:
````bash
git clone https://github.com/aganezov/RCK.git
````

then run the following command from the RCK source folder:
````bash
python setup.py install
```` 

**WARNING**: this installation does take care of python bindings for Gurobi. Please, refer to respective [subsection](#python-bindings-for-gurobi) on how that can be addressed.

### Gurobi 
[Gurobi](http://www.gurobi.com/) solver can be obtained from the official web site and installation procedure is also described there.
Gurobi requires a valid license to run. 
Licensing [information](http://www.gurobi.com/downloads/licenses/license-center) is provided on the official website, and is available for free for individual academic users.
More complicated setups with multi-user and cluster licenses are also available (and described on the official Gurobi website).
Contact your university IT support for more information about any complication with Gurobi licensing and setup.

RCK expects that Gurobi is installed on the machine in question.
RCK requires python bindings be installed (in the virtual environment, if you use it (which we recommend)). 
Refer to the next [subsection](#python-bindings-for-gurobi) for details on how this can be addressed.

##### Python bindings for Gurobi 
RCK requires python bindings be installed (in the virtual environment, if you use it (which we recommend)).
The following [documentation](https://www.gurobi.com/documentation/8.1/quickstart_windows/the_gurobi_python_interfac.html) of the Gurobi website explains how an installation of such bindings can be done.

Recommended way is via anaconda.
Regardless of whether conda is used for virtual environment, os just in general, the following command will install Python Guorbi bindings:
````bash
conda install -c gurobi gurobi
````

If not using conda, one needs to go to the Gurobi installation dir and locate the `setup.py` file and run the following command:
````bash
python setup.py install
````
Note that this way is deprecated by Gurobi.
 

### Executables
Installation of RCK adds several executables to your `PATH` (if using virtual environment, this executables will be accessible only when the environment is activated):
* `rck` - main executable that runs the RCK inference algorithm
* `rck-adj-x2rck` - conversion of SV prediction from several 3rd-party SV prediction tools (refer to respective [docs section](Adjacencies.md#converting-to-rck-format-from-sv-detection-tools) for more details)
* `rck-adj-process` - various processing options for RCK formatted adjacencies
* `rck-scnt-x2rck` - conversion of the clone- and allele-specific segment copy number predictions by 3rd-party tools (refer to respective [docs section](Segments.md#converting-to-rck-format-from-clone--and-allele-specific-inference-tools) for more details)
* `rck-scnt-process` - various processing options for RCK formatted segments, copy number, boundaries, etc. 