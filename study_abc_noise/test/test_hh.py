from study_abc_noise.model import HodgkinHuxleyModelVars as ModelVars
import matplotlib.pyplot as plt

# ModelVars.install_model()
mv = ModelVars()
gt = {"dc": 20, "membrane_dim": 10}
fig, axes = plt.subplots(nrows=4, sharex=True)
fig.set_size_inches((12, 8))
for _ in range(1):
    obs = mv.get_model()(gt)
    obs.plot(y="K", color="C1", ax=axes[0])
    obs.plot(y="Na", color="C0", ax=axes[1])
    obs = mv.get_model_noisy()(gt)
    obs.plot(y="K", color="C3", ax=axes[2], label="K noisy")
    obs.plot(y="Na", color="C2", ax=axes[3], label="Na noisy")
for ax in axes:
    ax.legend().set_visible(False)
axes[0].set_title("K")
axes[0].set_ylabel("K")
axes[1].set_title("Na")
axes[1].set_ylabel("Na")
plt.show()
