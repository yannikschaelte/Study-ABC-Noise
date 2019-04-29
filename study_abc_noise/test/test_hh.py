from study_abc_noise.model import HodgkinHuxleyModelVars as ModelVars
import matplotlib.pyplot as plt

# ModelVars.install_model()
mv = ModelVars(n_t=100)
gt = {"dc": 20, "membrane_dim": 10}
fig, axes = plt.subplots(nrows=2, sharex=True)
fig.set_size_inches((12, 8))
for _ in range(10):
    obs = mv.get_model()(gt)
    axes[0].plot(range(0, len(obs['K'])), obs["K"], 'x-')
    #obs.plot(y="Na", color="C0", ax=axes[1])
    obs = mv.add_noise(obs)
    # obs = mv.get_model_noisy()(gt)
    axes[1].plot(range(0, len(obs['K'])), obs['K'], 'x-')
    #obs.plot(y="K", color="C3", ax=axes[2], label="K noisy")
    #obs.plot(y="Na", color="C2", ax=axes[3], label="Na noisy")
for ax in axes:
    ax.legend().set_visible(False)
axes[0].set_title("K")
axes[0].set_ylabel("K")
#axes[1].set_title("Na")
#axes[1].set_ylabel("Na")
#plt.show()
plt.savefig("test/system_hh.png")
