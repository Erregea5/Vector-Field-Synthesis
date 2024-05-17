import matplotlib.pyplot as plt
import numpy as np

figures=[]

# num_colors=101
# colors=['']*num_colors
# for i in range(num_colors):
#   r=(i)/num_colors
#   colors[i]=(1,1-r,1-r)

def renderArrows(x,y,vx,vy,title,newFigure=True):
  scaling_factor =.05
  vx_scaled = vx * scaling_factor
  vy_scaled = vy * scaling_factor
  if newFigure:
    figures.append(plt.figure(figsize=(4, 4)))
  plt.quiver(x, y, vx_scaled, vy_scaled, pivot='middle', angles='xy', scale_units='xy', scale=1)
  plt.xlabel('X-axis')
  plt.ylabel('Y-axis')
  plt.title(title)
  plt.axis('equal')

def renderContours(x,y,z:np.ndarray,title):
  x=np.linspace(0,x.max(),z.shape[1])
  y=np.linspace(0,y.max(),z.shape[0])
  
  figures.append(plt.figure(figsize=(4, 4)))
  plt.contourf(x,y,z,levels=100,cmap='Reds')
  plt.colorbar()
  plt.xlabel('X-axis')
  plt.ylabel('Y-axis')
  plt.title(title)
  # plt.axis('equal')

def renderVertices(field,title,color=False,usePreviousFigure=False):
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
    renderContours(x,y,points,title+' color')
  else:
    renderArrows(x,y,vx,vy,title,not usePreviousFigure)

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

    renderContours(x,y,div.T,title+' divergence')
    renderContours(x,y,curl.T,title+' curl')
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

def renderField(field,title,render='vertices',usePreviousFigure=False):
  if render=='vertices':
    renderVertices(field,title,False,usePreviousFigure)
  elif render=='vertices-color':
    renderVertices(field,title,True,usePreviousFigure)
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