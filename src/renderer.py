import matplotlib.pyplot as plt
import numpy as np

figures=[]

def renderArrows(x,y,vx,vy,title):
  scaling_factor = 0.05
  vx_scaled = vx * scaling_factor
  vy_scaled = vy * scaling_factor

  figures.append(plt.figure(figsize=(4, 4)))
  plt.quiver(x, y, vx_scaled, vy_scaled, pivot='middle', angles='xy', scale_units='xy', scale=1)
  plt.xlabel('X-axis')
  plt.ylabel('Y-axis')
  plt.title(title)
  plt.axis('equal')

def renderContours(z:np.ndarray,title):
  x=np.linspace(0,1,z.shape[1])
  y=np.linspace(0,1,z.shape[0])
  figures.append(plt.figure(figsize=(4, 4)))
  plt.contourf(x,y,z,levels=100)
  plt.colorbar()
  plt.xlabel('X-axis')
  plt.ylabel('Y-axis')
  plt.title(title)
  plt.axis('equal')

def renderVertices(field,title,color=False):
  vertices = field.vertices
  x=np.zeros(len(vertices))
  y,vx,vy=x.copy(),x.copy(),x.copy()
  for i in range(len(vertices)):
    x[i]=vertices[i].pos[0]
    y[i]=vertices[i].pos[1]
    if vertices[i].dir is None:
      vx[i]=vy[i]=0
    else:
      vx[i]=vertices[i].dir[0]
      vy[i]=vertices[i].dir[1]
  if color:
    size_x=len(set(x))
    size_y=len(set(y))
    
    points=np.zeros((size_x+1,size_y+1))
    
    x_inv_step=size_x/(x.max()-x.min())
    y_inv_step=size_y/(y.max()-y.min())
    for i in range(len(x)):
      px=int((x[i]-x.min())*x_inv_step)
      py=int((y[i]-y.min())*y_inv_step)
      points[px,py]=(vx[i]**2+vy[i]**2)**.5
    renderContours(points,title+' color')
  else:
    renderArrows(x,y,vx,vy,title)

def renderJacobian(field,title,color=False):
  vertices = field.vertices
  x=np.zeros(len(vertices))
  y,dvx_dx,dvy_dx,dvx_dy,dvy_dy=x.copy(),x.copy(),x.copy(),x.copy(),x.copy()
  func=lambda x:x
  # if not color:
  #   func=vector.normed
  for i in range(len(vertices)):
    x[i]=vertices[i].pos[0]
    y[i]=vertices[i].pos[1]
    if vertices[i].jacobian is None:
      dvx_dx[i]=dvy_dx[i]=dvx_dy[i]=dvy_dy[i]=0
    else:
      dvx_dx[i],dvx_dy[i]=func(vertices[i].jacobian[0][0:2])
      dvy_dx[i],dvy_dy[i]=func(vertices[i].jacobian[1][0:2])

  if color:
    size_x=len(set(x))
    size_y=len(set(y))
    
    curl=np.zeros((size_x+1,size_y+1))
    div=curl.copy()
    x_inv_step=size_x/(x.max()-x.min())
    y_inv_step=size_y/(y.max()-y.min())
    for i in range(len(x)):
      px=int((x[i]-x.min())*x_inv_step)
      py=int((y[i]-y.min())*y_inv_step)
      div[px,py]=dvx_dx[i]+dvy_dy[i]
      curl[px,py]=dvy_dx[i]-dvx_dy[i]

    renderContours(div.T,title+' divergence')
    renderContours(curl.T,title+' curl')
  else:
    renderArrows(x,y,dvx_dx,dvx_dy,title+' dvx')
    renderArrows(x,y,dvy_dx,dvy_dy,title+' dvy')

def renderFaces(field,title):
  vertices = field.vertices
  x=np.zeros(len(vertices))
  y=x.copy()
  for i in range(len(vertices)):
    x[i]=vertices[i].pos[0]
    y[i]=vertices[i].pos[1]
    
  figures.append(plt.figure(figsize=(4, 4)))
  plt.triplot(x, y, field.faces)
  plt.plot(x, y, 'o')
  plt.xlabel('X-axis')
  plt.ylabel('Y-axis')
  plt.title(title)
  plt.axis('equal')

def renderField(field,title,render='vertices'):
  if render=='vertices':
    renderVertices(field,title)
  elif render=='vertices-color':
    renderVertices(field,title,True)
  elif render=='faces':
    renderFaces(field,title)
  elif render=='jacobian':
    renderJacobian(field,title)
  elif render=='jacobian-color':
    renderJacobian(field,title,True)

def renderLineChart(values:list[list],title,labels=None):
  figures.append(plt.figure(figsize=(4, 4)))
  if labels is None:
    for graph in values:
      plt.plot(graph)
  else:
    for graph,label in zip(values,labels):
      plt.plot(graph,label=label)
  plt.legend()
  plt.xlabel('X-axis')
  plt.ylabel('Y-axis')
  plt.title(title)
  # plt.axis('equal')

def renderHistogram(values,title,buckets=50):
  figures.append(plt.figure(figsize=(4, 4)))
  plt.hist(values,bins=buckets)
  plt.xlabel('X-axis')
  plt.ylabel('Y-axis')
  plt.title(title)

def close():
  for figure in figures:
    plt.close(figure)
  figure=[]

def show():
  plt.show(block=False)
  plt.pause(0.001) # Pause for interval seconds.
  input("hit[enter] to end.")
  close()