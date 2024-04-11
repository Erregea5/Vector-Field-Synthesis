from vector_field import vector,Field
import numpy as np
from renderer import renderField


def Nguyens_method(field):
  num_vertices=len(field.vertices)
  num_edges=[0]*num_vertices
  known_verts=[]


  for i in range(num_vertices):
    vec=field.vertices[i]
    if vec.jacobian and vec.dir:
      known_verts.append(i)
  
  for i in known_verts:
    vec=field.vertices[i]
    neighbors=field.edges[i]
    for j in neighbors:
      w_ij=field.edges[i][j]
      neighbor=field.vertices[j]
      aproximate_dir=vec.dir+vec.jacobian@(np.array(neighbor.pos)-np.array(vec.pos))
      if neighbor.dir!=None:
        for u in range(3):
          neighbor.dir[u]+=w_ij*aproximate_dir[u]
      else:
        neighbor.dir=[w_ij*x for x in aproximate_dir]
      num_edges[j]+=1

  for i in range(num_vertices):
    if num_edges[i]>0:
      field.vertices[i].dir=vector.normed(field.vertices[i].dir)
  renderField(field,'nguyen\'s intermediate')
  

def vectorSynthesis(field:Field,method='merged')->Field:
  '''field is a sparse field with or without jacobian information'''
  if len(field.edges)==0:
    field.calculate_edges()

  if method=='Nguyens':
    Nguyens_method(field)

  num_vertices=len(field.vertices)
  A=np.zeros((num_vertices*3,num_vertices*3))
  b=np.zeros(num_vertices*3)

  def coefficients_with_jacobian(i):
    '''equations: \n
        J_i*(sum(w_ij(p_j))-p_i) = sum(w_ij*v_j) - v_i                                         (jacobian at x) \n
        J_i*(sum(w_ij(p_j))-p_i) + sum(J_j*(p_j-p_i)) = sum(w_ij*v_j) + sum(v_j) - v_i*m       (jacobian at x and x+h_i)
    '''
    vec=field.vertices[i]
    edges=field.edges[i]
    p_i=np.array(vec.pos)
    h_i=-p_i
    val=np.zeros(3)
    m=1

    for j in edges:
      vec_j=field.vertices[j]
      w_ij=edges[j]
      p_j=np.array(vec_j.pos)
      h_i+=w_ij*p_j

      if vec_j.jacobian:
        val+=w_ij*(vec_j.jacobian@(p_j-p_i))
        w_ij+=w_ij
        m+=w_ij
      for u in range(3):
        A[i*3+u,j*3+u]=w_ij

    # if vec.jacobian:
    val+=vec.jacobian@h_i

    for u in range(3):
      A[i*3+u,i*3+u]=-m
      b[i*3+u]=val[u]

  def coefficients_with_neighbors(i):
    '''equations:\n
        sum(J_j*(p_j-p_i)) = sum(w_ij*v_j) - v_i*m                                  (jacobian at x+h_i)
        0 = sum( w_ij * v_j )  +  -1 * v_i
    '''
    vec=field.vertices[i]
    edges=field.edges[i]
    p_i=np.array(vec.pos)
    val=np.zeros(3)
    m=0

    for j in edges:
      vec_j=field.vertices[j]
      if vec_j.jacobian and vec_j.dir:#attempt at replicating nguyen's method in the linaer system
        w_ij=edges[j]
        p_j=np.array(vec_j.pos)      
        val+=w_ij*(vec_j.jacobian@(p_j-p_i))
        m+=w_ij
        for u in range(3):
          A[i*3+u,j*3+u]=w_ij
    
    if m==0:
      m=1
      for j in edges:
        w_ij=edges[j]
        for u in range(3):
          A[i*3+u,j*3+u]=w_ij

    for u in range(3):
      A[i*3+u,i*3+u]=-m
      b[i*3+u]=val[u]
  
  def coefficients_known(i):
    '''equation: v_i = 1 * v_i'''
    for u in range(3):
      A[i*3+u,i*3+u]=1
      b[i*3+u]=field.vertices[i].dir[u]

  for i in range(num_vertices):
    vec=field.vertices[i]
    if vec.dir:
      coefficients_known(i)
    elif vec.jacobian:
      coefficients_with_jacobian(i)
    else:
      coefficients_with_neighbors(i)

  x=np.linalg.solve(A,b)
  reconstructed_field=Field()
  reconstructed_field.faces=field.faces
  reconstructed_field.edges=field.edges
  for i in range(num_vertices):
    reconstructed_field.vertices.append(vector(field.vertices[i].pos))
    reconstructed_field.vertices[i].dir=vector.normed(x[i*3:(i+1)*3])

  return reconstructed_field
