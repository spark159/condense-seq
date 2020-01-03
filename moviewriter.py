"""
===========
MovieWriter
===========

This example uses a MovieWriter directly to grab individual frames and write
them to a file. This avoids any event loop integration, but has the advantage
of working with even the Agg backend. This is not recommended for use in an
interactive setting.

"""
# -*- noplot -*-

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import random

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=15, metadata=metadata)

fig = plt.figure()
#l, = plt.plot([], [], 'k-o')

#plt.xlim(-5, 5)
#plt.ylim(-5, 5)
#l, = plt.spy(np.zeros((2,2)))

x0, y0 = 0, 0

with writer.saving(fig, "writer_test.mp4", 100):
    for i in range(100):
        #x0 += 0.1 * np.random.randn()
        #y0 += 0.1 * np.random.randn()
        #l.set_data(x0, y0)
        A = np.zeros((2,2))
        for i in range(2):
            for j in range(2):
                value = random.choice([0,1])
                A[i][j] = value
        plt.spy(A)
        writer.grab_frame()
