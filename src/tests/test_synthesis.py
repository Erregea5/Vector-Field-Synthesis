import os,sys
src_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path=os.path.dirname(src_path)+'\\data\\'
sys.path.append(src_path)

from vector_field import vector,Field
from renderer import renderField,show
from synthesis import vectorSynthesis
import random

bnoise = Field(data_path+"bnoise.ply")
renderField(bnoise,'original')

sparse= Field()
sparse.edges=bnoise.edges
sparse.faces=bnoise.faces

for vertex in bnoise.vertices:
  sparse.vertices.append(vector(vertex.pos))
  if random.random()<.2:
    sparse.vertices[-1].dir=vertex.dir

renderField(sparse,'sparse')

reconstructed=vectorSynthesis(sparse)
renderField(reconstructed,'reconstructed')
show()