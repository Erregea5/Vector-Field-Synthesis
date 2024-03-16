import matplotlib.pyplot as plt
import numpy as np

def renderField(field,title):
  vertices = field.vertices

  x=np.zeros(len(vertices))
  y=x.copy()
  vx=x.copy()
  vy=x.copy()
  for i in range(len(vertices)):
    x[i]=vertices[i].pos[0]
    y[i]=vertices[i].pos[1]
    if vertices[i].dir is None:
      vx[i]=0
      vy[i]=0
    else:
      vx[i]=vertices[i].dir[0]
      vy[i]=vertices[i].dir[1]

  scaling_factor = 0.05
  vx_scaled = vx * scaling_factor
  vy_scaled = vy * scaling_factor

  plt.figure(figsize=(6, 6))
  plt.quiver(x, y, vx_scaled, vy_scaled, pivot='middle', angles='xy', scale_units='xy', scale=1)
  plt.xlabel('X-axis')
  plt.ylabel('Y-axis')
  plt.title(title)
  plt.axis('equal')

def show():
  plt.show()