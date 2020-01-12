# Study-ABC-Noise

This repository contains the code for the analysis on how to account for measurement noise and model error in sequential ABC. It supplements the paper "Efficient noise assessment in sequential ABC". TODO: Add reference.


# Installation


We require python3.6 or later (we used Anaconda3-2018.12), and an installation of the toolbox [pyabc](https://github.com/icb-dcm/pyabc). The used version 0.9.25 has been attached in the [pyabc](pyabc/) directory.

The repository is organized as a python package that can be easily installed by

```sh
pip3 install study_abc_noise
```

from the base folder. All requirements can be installed by

```sh
pip3 install -r requirements.txt
```

Using `pip freeze`, we have here logged all versions of the packages we used, though also more recent versions may work.


## Organization

 
Globally, we have attached the code of the package [pyabc](https://github.com/icb-dcm/pyabc) in the version 0.9.25 that we used.


# Running the code


Many of the code blocks should run out of the box. For some, however absolute paths were necessary, which will need to be adapted.


## Figure guide


Here's where to find the code responsible for the figures contained in the paper and supplement.

Main manuscript:

* Motivation figure for wrong noise assessment for conversion reaction model
