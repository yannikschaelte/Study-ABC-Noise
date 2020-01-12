# Study-ABC-Noise

This repository contains the code for the analysis on how to account for measurement noise and model error in sequential ABC. It supplements the paper "Efficient noise assessment in sequential ABC". TODO: Add reference.


# Installation


We require python3.6 or later (we used Anaconda3-2018.12), and an installation of the toolbox [pyabc](https://github.com/icb-dcm/pyabc). The used version 0.9.25 has been attached in the [pyabc](pyabc/) directory. This is where most of the code relevant for this work has gone into.

The repository is organized as a python package that can be easily installed by

```sh
pip3 install study_abc_noise
```

from the base folder. All requirements can be installed by

```sh
pip3 install -r requirements.txt
```

We have here logged all versions of the packages we used, using `pip freeze`, though also more recent versions may work.


## Organization



### Folder structure


In the base directory, there are notably the folders `pyabc` and `study_abc_noise`. The former just contains a copy of the pyabc package in the version used, the latter we will describe now in more detail.

* The [model](study_abc_noise/model) folder contains the models that have been used.
* The file [vars](study_abc_noise/vars.py) contains some convenience code that allows handling all models in a similar fashion.
* The `application_...` folders contain application projects
 

### In each folder


* In each folder, there is a `README.md` file containing further information.
* Usually, in each folder there is a `script_run.py` or `run.ipynb` or similar file, which contains the code for executing the simulation. Further, there are `script_visualization.py` or `viz.ipynb` or similar files containing code for creating figures and other analysis.
* As far as possible, we have made use of jupyter notebooks, which allow to easily combine code and nalysis.


# Running the code


Many of the code blocks should run out of the box. For some, however absolute paths were necessary, which will need to be adapted.

For the sampling process we used two samplers: A `pyabc.sampler.MulticoreEvalParallelSampler` and a `pyabc.sampler.RedisEvalParallelSampler`. The former should run out of the box and can exploit parallelism on a single machine. The latter however, which can exploit cluster infrastructure via a redis server, needs a specific setup depending on the cluster infrastructure. See [here](https://pyabc.rtfd.io/en/latest/sampler.html) for further information.


## Figure guide


Here's where to find the code responsible for the figures contained in the paper and supplement.

Main manuscript:

* Motivation figure for wrong noise assessment for conversion reaction model
