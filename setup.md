# Installation of Conda

Conda is an open-source tool that provides package, dependency, and environment management for any programming language. To install conda, we must first pick the right installer for us. Below we will demonstrate how to install **Anaconda Distribution**, a full featured installer of Conda.


1. [Download the installer by choosing the proper installer for your machine.](https://www.anaconda.com/download/)
2. [Verify your installer hashes using SHA-256 checksums.](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#hash-verification)
3. Install the installer:
	- Windows: Double-click the `.exe` file. Follow the instructions on the screen. For a detailed reference, please read [this page](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html#installing-on-windows).
	- macOS: double-click the `.pkg` file. Follow the instructions on the screen.For a detailed reference, please read [this page](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html#installing-on-macos).
	- Linux: In your terminal window, run: `bash Anaconda-latest-Linux-x86_64.sh`. Follow the prompts on the installer screens. For a detailed reference, please read [this page](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html#installing-on-linux).


# Managing Environments

With conda, you can create, export, list, and update environments that have different versions of Python and/or packages installed in them. Below we will demonstrate how to create an environment for this tutorial on macOS/Linux. Use the **terminal** for the following steps. For a detailed reference, please read [this page](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

1. To create an environment with the latest version of python: `conda create --name <my-env>`. Replace `<my-env>` with the name of your environment.
2. When conda asks you to proceed, type `y`.

# Managing Packages

Before installing packages, use the **terminal** for the following steps:

1. Activate the environment you just created: `conda activate <my-env>`
2. (Optional) Check the Python version: `python --version`
3. (Optional) Check the pip version: ` pip --version`

The installation steps below are all run in the **terminal**.

## Installing Jupyter Notebook
The Jupyter Notebook is a web application for computational documents so that our code can produce rich, interactive output. 

1. To install notebook: `pip install notebook`
2. To run notebook: `jupyter notebook`

## Installing packages for analyzing genomics data

`pip install anndata` ([ocumentatoin of anndata](https://anndata.readthedocs.io/en/latest/index.html))

`pip install scanpy` ([Documentatoin of Scanpy](https://scanpy.readthedocs.io/en/stable/))


## Installing packages for computation with genomics data

### STAN

`pip install matplotlib`

`pip install numpy`

`pip install pandas`

`pip install Pillow`

`pip install seaborn`

`pip install scipy`

`pip install sklearn`

We will compare STAN with [decoupler](https://decoupler-py.readthedocs.io/en/latest/index.html).

`pip install decoupler`

### xxx (Parham)

### xxx (Xiaojun)

