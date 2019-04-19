from study_abc_noise.model import HodgkinHuxleyModelVars as ModelVars
import matplotlib.pyplot as plt


mv = ModelVars()
gt = {"dc": 20, "membrane_dim": 10}
fig, axes = plt.subplots(nrows=2, sharex=True)
fig.set_size_inches((12, 8))
for _ in range(10):
    obs = mv.get_model()(gt)
    obs.plot(y="K", color="C1", ax=axes[0])
    obs.plot(y="Na", color="C0", ax=axes[1])
for ax in axes:
    ax.legend().set_visible(False)
axes[0].set_title("K")
axes[0].set_ylabel("K")
axes[1].set_title("Na")
axes[1].set_ylabel("Na")
plt.show()
