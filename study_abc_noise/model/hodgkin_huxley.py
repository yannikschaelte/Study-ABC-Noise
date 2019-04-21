from ..vars import ModelVars
import requests
import os
from zipfile import ZipFile
import subprocess
from io import BytesIO
import numpy as np
import pandas as pd


model_folder = "hodgkin_huxley"
executable = os.path.join(model_folder, "ModelDBFolder", "HH_run")


class HodgkinHuxleyModelVars(ModelVars):

    def __init__(self):
        self.noise_std = 0.005

    def get_model(self):
        def model(p):
            return simulate(**p)
        return model

    def get_model_noisy(self):
        model = self.get_model()
        def model_noisy(p):
            ret = model(p)
            ret['Na'] += self.noise_std * np.random.randn(10001)
            ret['K'] += self.noise_std * np.random.randn(10001)
            return ret
        return model_noisy

    @staticmethod
    def install_model():
        url = ("https://senselab.med.yale.edu/modeldb/"
               "eavBinDown.cshtml?o=128502&a=23&mime=application/zip")
        if not os.path.isdir(model_folder):
            os.mkdir(model_folder)
        req = requests.request("GET", url)
        archive = ZipFile(BytesIO(req.content))
        archive.extractall(model_folder)
        ret = subprocess.run(
            ['make', 'HH_run'],
            cwd=os.path.join(model_folder, "ModelDBFolder"))
        print(f"The executable location is {executable}.")


def simulate(model=2, membrane_dim=10, time_steps=1e4, time_step_size=0.01,
             isi=100, dc=20, noise=0, sine_amplitude=0, sine_frequency=0,
             voltage_clamp=0, data_to_print=1, rng_seed=None):
    """
    Simulate the SDE ion channel model defined in a nexternal fortran
    simulator.

    Returns
    -------

    pandas.DataFrame
        index: t, time
        columns: V, Na, K
        V: Voltage
        Na, K: Proportion of open channels
    """
    if rng_seed is None:
        rng_seed = np.random.randint(2**32 - 2) + 1
    membrane_area = membrane_dim**2
    # print("begin sim")
    re = subprocess.run(
        [executable, str(model),
         # long floats cannot be handled
         f"{membrane_area:.5f}", str(time_steps),
         str(time_step_size), str(isi), f"{dc:.5f}",
         str(noise), str(sine_amplitude), str(sine_frequency),
         str(voltage_clamp), str(data_to_print),
         str(rng_seed)],
        stdout=subprocess.PIPE)
    df = pd.read_csv(BytesIO(re.stdout),
                     delim_whitespace=True,
                     header=None, index_col=0,
                     names=["t", "V", "Na", "K"])
    # print("end sim")
    return df
