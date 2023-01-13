import matplotlib.pyplot as plt
import numpy as np

balls = np.loadtxt('build/balls.txt')

cx = balls[:,0]
cy = balls[:,1]
px = balls[:,3]
py = balls[:,4]

plt.scatter(px, py, marker='.', c='k', label='surface points')
plt.scatter(cx, cy, marker='.', c='tab:red', label='medial ball centers')
ax = plt.gca()
for b in balls:
    ax.add_patch(plt.Circle((b[0], b[1]), b[6], edgecolor='tab:red', facecolor='none', alpha=0.2))
plt.axis('image')
plt.xlim([-1.5, +2.0])
plt.ylim([-1.5, +1.4])
plt.legend()
plt.show()
