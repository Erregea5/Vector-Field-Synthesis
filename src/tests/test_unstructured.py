import os,sys
src_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path=os.path.dirname(src_path)+'\\data\\'
sys.path.append(src_path)

from vector_field import vector,Field
from renderer import renderField,show
from synthesis import vectorSynthesis
import random

original = Field(data_path+"bnoise.ply")

expected=Field()
percent_vert_present=.15
for vertex in original.vertices:
  if random.random()<percent_vert_present:
    expected.vertices.append(vertex.copy())

renderField(expected,'original')
renderField(expected,'edges',render='faces')

sparse= Field()
percent_present=.3
for vertex in expected.vertices:
  sparse.vertices.append(vector(vertex.pos))
  if random.random()<percent_present:
    sparse.vertices[-1].dir=vertex.dir

renderField(sparse,'sparse')

reconstructed=vectorSynthesis(sparse)
renderField(reconstructed,'reconstructed')

error_field,err=Field.get_Error(reconstructed,expected)
print('error: ',err)
renderField(error_field,'error')

#delaunay triangulation
#try to calculate jacobian
#test the two methods
#try to make unstructured gird
#combine the two
#2 probing questions - 
# When working with an unstructured grid what are your thoughts on treating it like a sparse field in a rectilinear grid for simplified calculation?
# What are some examples in the real worrld where we'd deal with an unstructured grid?
# What function should i use to calculate error, and why?
#line integral convolution 

show()