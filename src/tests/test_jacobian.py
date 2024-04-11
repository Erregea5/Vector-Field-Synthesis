import os,sys
src_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path=os.path.dirname(src_path)+'\\data\\'
sys.path.append(src_path)

from vector_field import vector,Field
from renderer import renderField,show

original = Field(data_path+"bnoise.ply")
expected = original
expected.calculate_jacobian()

renderField(expected,'original')
renderField(expected,'jacobian arrows',render='jacobian')
renderField(expected,'jacobian colors',render='jacobian-color')

show()

#create elevator pitch for project
#degub!!
#set seed for sparse field creation and define point radius for further apart points
