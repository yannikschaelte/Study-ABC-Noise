from ..vars import ModelVars
import requests
import os
from zipfile import ZipFile
import subprocess
from io import BytesIO
import numpy as np
import pandas as pd
import copy
import pyabc


model_folder = ".hodgkin_huxley"
executable = os.path.join(model_folder, "ModelDBFolder", "HH_run")


class HodgkinHuxleyModelVars(ModelVars):

    def __init__(self, p_true: dict = None, n_acc: int = 100, n_t: int = 50,
                 noise_std: float = 0.05, basedir='.'):
        if p_true is None:
            p_true = {'dc': 20, 'membrane_dim': 10}
        super().__init__(p_true = p_true, n_acc=n_acc)
        self.noise_std = noise_std
        self.limits = {'dc': (2, 30), 'membrane_dim': (1, 12)}
        self.time_steps = 10001
        self.n_t = n_t
        self.executable = os.path.join(basedir, model_folder, "ModelDBFolder", "HH_run")
    
    def get_obs_times(self):
        return [int(i) for i in np.linspace(0, self.time_steps - 1, self.n_t)]

    def get_id(self):
        return f"hodgkin_huxley_{self.noise_std}"

    def get_prior(self):
        return pyabc.Distribution(
            **{key: pyabc.RV('uniform', bounds[0], bounds[1])
               for key, bounds in self.limits.items()})

    def get_distance(self):
        def l2(x, y):
            return np.sum(np.power( (x['K'] - y['K']), 2))
        return l2

    def get_kernel(self):
        kernel = pyabc.distance.IndependentNormalKernel(
            var=self.noise_std**2 * np.ones(self.n_t),
            pdf_max=self.pdf_max,
            keys=['K'])
        return kernel

    def get_model(self):
        def model(p):
            val = simulate(**p, executable=self.executable)
            ret = {'K': val.reset_index()['K'][self.get_obs_times()]}
            return ret
        return model

    def get_model_noisy(self):
        model = self.get_model()
        def model_noisy(p):
            ret = model(p)
            ret = self.add_noise(ret)
            return ret
        return model_noisy
    
    def add_noise_to_data(self, y):
        noisy_y = copy.deepcopy(y)
        noisy_y['K'] += self.noise_std * np.random.randn(self.n_t)
        return noisy_y

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
             voltage_clamp=0, data_to_print=1, rng_seed=None, executable=executable):
    """
    Simulate the SDE ion channel model defined in an external fortran
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
