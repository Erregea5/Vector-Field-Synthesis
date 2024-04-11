import os,sys
src_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path=os.path.dirname(src_path)+'\\data\\'
sys.path.append(src_path)

from vector_field import vector,Field
from renderer import renderField,show
from synthesis import vectorSynthesis
import random

original = Field(data_path+"bnoise.ply")
expected = original
expected=Field()
percent_vert_present=.15
for vertex in original.vertices:
  if random.random()<percent_vert_present:
    expected.vertices.append(vertex.copy())

expected.calculate_jacobian()
renderField(expected,'original')

def copy():
  field=Field()
  field.edges=expected.edges
  field.faces=expected.faces
  return field
percent_present=.05
sparse=copy()
sparse_j=copy()

for vertex in expected.vertices:
  sparse.vertices.append(vector(vertex.pos))
  sparse_j.vertices.append(vector(vertex.pos))
  sparse_j.vertices[-1].jacobian=vertex.jacobian
  if random.random()<percent_present:
    sparse.vertices[-1].dir=vertex.dir
    sparse_j.vertices[-1].dir=vertex.dir

renderField(sparse_j,'sparse')

reconstructed=vectorSynthesis(sparse)
reconstructed_j=vectorSynthesis(sparse_j)
renderField(reconstructed,'reconstructed without jacobian')
renderField(reconstructed_j,'reconstructed with jacobian')

error_field,err=Field.get_Error(reconstructed,expected)
error_field_j,err_j=Field.get_Error(reconstructed_j,expected)
dif_field,dif=Field.get_Error(reconstructed_j,reconstructed)
print('error: ',err)
renderField(error_field,'regular error')
print('error with jacobian: ',err_j)
renderField(error_field_j,'error with jacobian')
print('difference: ',dif)
renderField(dif_field,'difference')

show()