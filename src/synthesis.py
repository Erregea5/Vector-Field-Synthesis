from vector_field import vector,Field
import numpy as np
from renderer import renderField,close,show


def Nguyens_method(field,normalized):
  num_vertices=len(field.vertices)
  num_edges=[0]*num_vertices
  known_verts=set()


  for i in range(num_vertices):
    vec=field.vertices[i]
    if vec.jacobian is not None and vec.dir:
      known_verts.add(i)
  
  for i in known_verts:
    vec=field.vertices[i]
    neighbors=field.edges[i]
    for j in neighbors:
      if j in known_verts:
        continue
      w_ij=field.edges[i][j]
      neighbor=field.vertices[j]
      aproximate_dir=vec.dir+vec.jacobian@(np.array(neighbor.pos)-np.array(vec.pos))
      if neighbor.dir is not None:
        # for u in range(3):
        neighbor.dir+=w_ij*aproximate_dir
      else:
        neighbor.dir=w_ij*aproximate_dir
      num_edges[j]+=w_ij

  for i in range(num_vertices):
    vec=field.vertices[i]
    vec.jacobian=None
    if normalized:
      vec.dir=vector.normed(vec.dir)
    elif num_edges[i]>0:
      vec.dir/=num_edges[i]
      vec.dir=[x for x in vec.dir]
  
def vectorSynthesisWithAugmentedMatrix(field:Field,method='curl',normalized=False)->Field:
  if len(field.edges)==0:
    field.calculate_edges()

  num_physics=0
  for i in range(len(field.vertices)):
    if field.vertices[i].jacobian is not None: 
      num_physics+=1
  if method=='jacobian':
    num_physics*=2

  num_vertices=len(field.vertices)
  A=np.zeros((num_vertices*3+num_physics,num_vertices*3))
  b=np.zeros(num_vertices*3+num_physics)

  physics_cnt=0

  def set_curl(eq,inv_dP,curl,i,j):
    w_ij=field.edges[i][j]
    inv_dP*=w_ij
    eq[i*3+1]+=inv_dP[0]
    eq[i*3]+=-inv_dP[1]
    eq[j*3+1]+=-inv_dP[0]
    eq[j*3]+=inv_dP[1]

    b[num_vertices*3+physics_cnt]+=curl*w_ij
  
  def set_div(eq,inv_dP,div,i,j):
    w_ij=field.edges[i][j]
    inv_dP*=w_ij
    eq[i*3]+=inv_dP[0]
    eq[i*3+1]+=inv_dP[1]
    eq[j*3]+=-inv_dP[0]
    eq[j*3+1]+=-inv_dP[1]

    b[num_vertices*3+physics_cnt]+=div*w_ij

  def set_jacobian(i,j):
    nonlocal physics_cnt
    vec=field.vertices[i]
    vec_j=field.vertices[j]
    p_i=np.array(vec.pos)
    p_j=np.array(vec_j.pos)
    val=vec.jacobian@(p_j-p_i)
    w_ij=field.edges[i][j]

    for u in range(2):
      k=num_vertices*3+physics_cnt+u
      A[k,i*3+u]+=-w_ij
      A[k,j*3+u]+=w_ij
      b[k]+=val[u]*w_ij

  def coefficients_with_physics(i):
    nonlocal physics_cnt
    vec=field.vertices[i]
    edges=field.edges[i]
    p_i=np.array(vec.pos)
    
    curl=vec.jacobian[1,0]-vec.jacobian[0,1]
    div=vec.jacobian[0,0]+vec.jacobian[1,1]
    for j in edges:
      vec_j=field.vertices[j]
      p_j=np.array(vec_j.pos)
      dP=p_i-p_j
      inv_dP=dP
      for u in range(len(dP)):
        inv_dP[u]=1/dP[u] if dP[u]!=0 else 0
      
      eq=A[num_vertices*3+physics_cnt]
      if method=='curl':
        set_curl(eq,inv_dP,curl,i,j)
      elif method=='div':
        set_div(eq,inv_dP,div,i,j)
      else: set_jacobian(i,j)
    if method=='jacobian':
      physics_cnt+=2
    else:
      physics_cnt+=1

  def coefficients_with_neighbors(i):
    '''equations:\n
        sum(J_j*(p_j-p_i)) = sum(w_ij*v_j) - v_i*m                                  (jacobian at x+h_i)
        0 = sum( w_ij * v_j )  +  -1 * v_i
    '''
    edges=field.edges[i]
    for j in edges:
      w_ij=edges[j]
      for u in range(3):
        A[i*3+u,j*3+u]+=w_ij

    for u in range(3):
      A[i*3+u,i*3+u]=-1
  
  def coefficients_known(i):
    '''equation: v_i = 1 * v_i'''
    for u in range(3):
      A[i*3+u,i*3+u]=1
      b[i*3+u]=field.vertices[i].dir[u]
  
  for i in range(num_vertices):
    vec=field.vertices[i]
    if vec.dir:
      coefficients_known(i)
    else:
      coefficients_with_neighbors(i)
    if vec.jacobian is not None:
      coefficients_with_physics(i)
    

  x=np.linalg.lstsq(A,b,rcond=None)[0]
  reconstructed_field=Field()
  reconstructed_field.faces=field.faces
  reconstructed_field.edges=field.edges
  for i in range(num_vertices):
    reconstructed_field.vertices.append(field.vertices[i].copy())
    dir=[c for c in x[i*3:(i+1)*3]]
    if normalized:
      dir=vector.normed(dir)
    reconstructed_field.vertices[i].dir=dir
    reconstructed_field.vertices[i].jacobian=field.vertices[i].jacobian

  return reconstructed_field

def vectorSynthesis(field:Field,method='merged',actual=None,normalized=False)->Field:
  '''field is a sparse field with or without jacobian information'''
  if len(field.edges)==0:
    field.calculate_edges()

  if method=='Nguyens':
    Nguyens_method(field,normalized)

  num_vertices=len(field.vertices)
  A=np.zeros((num_vertices*3,num_vertices*3))
  b=np.zeros(num_vertices*3)

  def coefficients_with_jacobian(i):
    '''equations: \n
        J_i*(sum(w_ij(p_j))-p_i) = sum(w_ij*v_j) - v_i                                         (jacobian at x) \n
    '''
    vec=field.vertices[i]
    edges=field.edges[i]
    p_i=np.array(vec.pos)
    h_i=-p_i

    for j in edges:
      vec_j=field.vertices[j]
      w_ij=edges[j]
      p_j=np.array(vec_j.pos)
      h_i+=p_j*w_ij
      for u in range(3):
        A[i*3+u,j*3+u]+=w_ij

    val=vec.jacobian@h_i
    for u in range(3):
      A[i*3+u,i*3+u]+=-1
      b[i*3+u]+=val[u]

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
      if vec_j.jacobian is not None and vec_j.dir:
        w_ij=edges[j]
        p_j=np.array(vec_j.pos)      
        val+=w_ij*(vec_j.jacobian@(p_j-p_i))
        m+=w_ij
        for u in range(3):
          A[i*3+u,j*3+u]+=w_ij
    
    if m==0:
      m=1
      for j in edges:
        w_ij=edges[j]
        for u in range(3):
          A[i*3+u,j*3+u]+=w_ij

    for u in range(3):
      A[i*3+u,i*3+u]=-m
      b[i*3+u]=val[u]
  
  def coefficients_known(i):
    '''equation: v_i = 1 * v_i'''
    for u in range(3):
      A[i*3+u,i*3+u]=1
      b[i*3+u]=field.vertices[i].dir[u]

  def check(i,reconstructed):
    vec_i=reconstructed.vertices[i]
    if vec_i.jacobian is not None:
      real_vec_i=actual.vertices[i]
      jj=0
      for j in reconstructed.edges[i]:
        jj=j; break
      vec_j=reconstructed.vertices[jj]
      real_vec_j=actual.vertices[jj]
      h_i=np.array(vec_i.pos)-vec_j.pos
      dif= sum((vec_i.jacobian@h_i - (np.array(vec_i.dir)-vec_j.dir))**2)
      if dif>.1:
        print(dif)

  for i in range(num_vertices):
    vec=field.vertices[i]
    if vec.dir:
      coefficients_known(i)
    elif vec.jacobian is not None:
      coefficients_with_neighbors(i)
      coefficients_with_jacobian(i)
    else:
      coefficients_with_neighbors(i)

  err=False
  try:
    x=np.linalg.solve(A,b)
  except np.linalg.LinAlgError:
    print('Error: Singular Matrix')
    x=np.linalg.lstsq(A,b)[0]
    err=True
  reconstructed_field=Field()
  reconstructed_field.faces=field.faces
  reconstructed_field.edges=field.edges
  for i in range(num_vertices):
    reconstructed_field.vertices.append(field.vertices[i].copy())
    dir=[c for c in x[i*3:(i+1)*3]]
    if normalized:
      dir=vector.normed(dir)
    reconstructed_field.vertices[i].dir=dir
    reconstructed_field.vertices[i].jacobian=field.vertices[i].jacobian

  # for i in range(num_vertices):
  #   check(i,reconstructed_field)

  if err:
    renderField(reconstructed_field,'fail')
    show()
  return reconstructed_field
