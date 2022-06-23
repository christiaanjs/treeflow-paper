import pickle

with open("vi-trace.pickle", "rb") as f:
    res = pickle.load(f)

import matplotlib.pyplot as plt

plt.plot(res.loss.numpy())
plt.show()
