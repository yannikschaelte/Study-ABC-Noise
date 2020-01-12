# Application to an SDE model of Hodgkin-Huxley neuron spiking

The model is taken from:

Goldwyn, Joshua H., Nikita S. Imennov, Michael Famulare, and Eric Shea-Brown. “Stochastic Differential Equation Models for Ion Channel Noise in Hodgkin-Huxley Neurons.” Physical Review E 83, no. 4 (2011): 041908. doi:10.1103/PhysRevE.83.041908.


## Notes

The model is written in fortran and must be compiled. The compiled files are in the `.hodgkin_huxley` folder, however it may be necessary to recompile. To do so, we have provided a convenience function handling the download and installation. Just open a python shell and execute:

```python
from study_abc_noise.model import HodgkinHuxleyModelVars as MV
MV.install_model()
```

It will install to the working directory, paths may need to be adapted thereafter.
