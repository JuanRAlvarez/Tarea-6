"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
datos = np.loadtxt('3cuerpos.dat')
x1 = datos[:,0]
y1 = datos[:,1]
x2 = datos[:,3]
y2 = datos[:,4]
x3 = datos[:,6]
y3 = datos[:,7]
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(x1[np.argmin(x1)], x1[np.argmax(x1)]), ylim=(y1[np.argmin(y1)], y1[np.argmax(y1)]))
line, = ax.plot([], [], 'or')

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(j):
    i = 20*j
    x = [x1[i],x2[i],x3[i]]
    y = [y1[i],y2[i],y3[i]]
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(x1)/10, interval=1, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
plt.show()
