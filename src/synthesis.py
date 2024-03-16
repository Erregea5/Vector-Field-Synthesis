from vector_field import vector,Field
import numpy as np

def weight(field:Field,i,j):
  if j in field.edges[i]:
    return 1.0/field.edges[i][j]
  return 0

def normed(arr):
  s=sum([(x)**2 for x in arr])**.5
  if s==0:s=1
  return [x/s for x in arr]

# vector field has 
def vectorSynthesis(sparse_field:Field)->Field:
  if len(sparse_field.edges)==0:
    sparse_field.calculate_edges()

  num_vertices=len(sparse_field.vertices)
  A=np.zeros((num_vertices*3,num_vertices*3))
  b=np.zeros(num_vertices*3)

  for i in range(num_vertices):
    dir=sparse_field.vertices[i].dir

    # if vector is in sparse field
    if dir is not None:
      for u in range(3):
        A[i*3+u,i*3+u]=1
        b[i*3+u]=dir[u]
      continue;
    else: 
      for u in range(3):
        A[i*3+u,i*3+u]=-1

    if i+1>num_vertices:
      break
    for j in range(i+1,num_vertices):
      w_ij=weight(sparse_field,i,j)
      for u in range(3):
        A[i*3+u,j*3+u]=w_ij

      dir=sparse_field.vertices[j].dir
      if dir is None:
        for u in range(3):
          A[j*3+u,i*3+u]=w_ij
  
  for i in range(num_vertices*3):
    if A[i,i]==1:
      continue
    A[i]=A[i]/(A[i].sum()+1)
    A[i,i]=-1

  x=np.linalg.solve(A,b)
  approximate_field=Field()
  approximate_field.faces=sparse_field.faces
  approximate_field.edges=sparse_field.edges
  for i in range(num_vertices):
    approximate_field.vertices.append(vector(sparse_field.vertices[i].pos))
    approximate_field.vertices[i].dir=([x[i*3+u] for u in range(3)])

  return approximate_field