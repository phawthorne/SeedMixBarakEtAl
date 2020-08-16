import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('random_samples/means.csv')
xcol = "SpecsInMix"
ycols = ["PD", "Shannon", "Cost","ConsVal","Bloom"]

fig, axes = plt.subplots(3, 2, figsize=(8, 10))

fig.suptitle("Means")
for i, ycol in enumerate(ycols):
    r = int(np.floor(i/2))
    c = int(i % 2)
    print(f'{r} {c}')
    ax = axes[r][c]
    df.plot.scatter(xcol, ycol, ax=ax)

axes[2][1].axis('off')

plt.savefig('random_samples/means.pdf')
plt.close()



df = pd.read_csv('random_samples/stddevs.csv')
xcol = "SpecsInMix"
ycols = ["PD", "Shannon", "Cost","ConsVal","Bloom"]

fig, axes = plt.subplots(3, 2, figsize=(8, 10))

fig.suptitle("Std Devs")
for i, ycol in enumerate(ycols):
    r = int(np.floor(i/2))
    c = int(i % 2)
    print(f'{r} {c}')
    ax = axes[r][c]
    df.plot.scatter(xcol, ycol, ax=ax)

axes[2][1].axis('off')

plt.savefig('random_samples/stddevs.pdf')
plt.close()


df = pd.read_csv('random_samples/medians.csv')
xcol = "SpecsInMix"
ycols = ["PD", "Shannon", "Cost","ConsVal","Bloom"]

fig, axes = plt.subplots(3, 2, figsize=(8, 10))

fig.suptitle("Medians")
for i, ycol in enumerate(ycols):
    r = int(np.floor(i/2))
    c = int(i % 2)
    print(f'{r} {c}')
    ax = axes[r][c]
    df.plot.scatter(xcol, ycol, ax=ax)

axes[2][1].axis('off')

plt.savefig('random_samples/medians.pdf')
plt.close()
