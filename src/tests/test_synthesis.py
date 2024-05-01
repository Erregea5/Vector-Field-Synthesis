import os,sys
src_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path=os.path.dirname(src_path)+'\\data\\'
sys.path.append(src_path)

from vector_field import vector,Field
from renderer import renderField,show
from synthesis import vectorSynthesis
import random

expected = Field(data_path+"bnoise.ply")
renderField(expected,'original')

percent_present=.05
sparse= Field()
sparse.edges=expected.edges
sparse.faces=expected.faces

for vertex in expected.vertices:
  sparse.vertices.append(vector(vertex.pos))
  if random.random()<percent_present:
    sparse.vertices[-1].dir=vertex.dir

renderField(sparse,'sparse')

reconstructed=vectorSynthesis(sparse)
renderField(reconstructed,'reconstructed')

error_field,err=Field.get_direction_error(reconstructed,expected)
print('error: ',err)
renderField(error_field,'error')

show()