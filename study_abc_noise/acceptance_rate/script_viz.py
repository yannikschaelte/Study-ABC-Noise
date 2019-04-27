import cloudpickle as pickle
import matplotlib.pyplot as plt


with open('dim.dat', 'rb') as f:
    val = pickle.load(f)

list_pdf_max = val[0]
list_pdf_max_min = val[1]
list_acc_rate_prior = val[2]
list_acc_rate_prior_min = val[3]
list_acc_rate_posterior = val[4]
list_acc_rate_posterior_min = val[5]
list_n_t = [4, 7, 10, 13, 16, 19, 22, 25, 28]
print(val)

_, ax = plt.subplots()
ax.semilogy(list_n_t, list_pdf_max, 'x-', label="$pdf_{max}$")
ax.semilogy(list_n_t, list_pdf_max_min, 'x-', label="$pdf_{opt}$")
ax.set_xlabel("#Data")
ax.set_ylabel("Maximum density")
ax.legend()
plt.savefig("pdf_maxs.png")


_, ax = plt.subplots()
ax.semilogy(list_n_t, list_acc_rate_prior, 'x-', label="prior")
ax.semilogy(list_n_t, list_acc_rate_prior_min, 'x-', label="prior opt")
ax.semilogy(list_n_t, list_acc_rate_posterior, 'x-', label="posterior")
ax.semilogy(list_n_t, list_acc_rate_posterior_min, 'x-', label="posterior opt")
ax.set_xlabel("#Data")
ax.set_ylabel("Acceptance rate")
ax.legend()
plt.savefig("acceptance_rates.png")

