import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Define the sigmoid function
def sigmoid(x, k):
    return (1 / (1 + np.exp(-k * x))) - 0.5

# Define the equation to find the intersection
def intersection_eq(x):
    return 2 * sigmoid(x, 1) - 0.5

x_intersect = fsolve(intersection_eq, 0)  # initial guess is 0
y_intersect = 0.5

x = np.linspace(0, 7, 400)  # Only positive x-values

# Sigmoid with k=1 (less steep)
y1 = 2 * sigmoid(x, 1)

# Sigmoid with k=2 (more steep)
y2 = 2.5 * sigmoid(x, 2)

# Horizontal straight line at y=0.5
y_horizontal = 0.5 * np.ones_like(x)

# Fourth line that follows the minimum of the horizontal straight line and k=1 sigmoid
y_min = np.minimum(y1, y_horizontal)

fig, ax = plt.subplots()
ax.plot(x, y1, label='J', color = 'black')
ax.plot(x, y2, label='Jtot', color = 'black')
ax.plot(x, y_horizontal, label='Jv', color = 'black', linestyle='--')
ax.plot(x, y_min, label='Jco2+o2', linestyle='-.', linewidth=1.5, color ='green')
ax.scatter(x_intersect, y_intersect, color='red', zorder=5, label='Coordination hypothesis')

ax.fill_between(x, y2, y_min, where=(y2 > y_min), color='pink', alpha=0.5)
ax.fill_between(x, y1, y_min, where=(y1 > y_min), facecolor='none', edgecolor='maroon', hatch='//', alpha=0.7)

ax.axhline(0, color='black', linewidth=0.5)
ax.axvline(0, color='black', linewidth=0.5)
ax.set_title("Sigmoid Functions and Horizontal Line with Intersection")
ax.legend()
#plt.grid(True, which='both', linestyle='--', linewidth=0.5)
ax.set_xlim(0, 4)
ax.set_ylim(0, 2)
ax.set_xlabel("PPFD (mol m2 s1)")
ax.set_ylabel("Electron Flux (!mol m!2 s!1)")
ax.set_xticks([])
ax.set_yticks([])

plt.show()





