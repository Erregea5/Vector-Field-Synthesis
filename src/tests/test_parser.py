import os,sys
src_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path=os.path.dirname(src_path)+'\\data\\'
sys.path.append(src_path)

from vector_field import Field
from renderer import renderField,show

bnoise = Field(data_path+"bnoise.ply")

renderField(bnoise,'field')
renderField(bnoise,'original faces',render='faces')
bnoise.faces=[]
bnoise.calculate_faces()
renderField(bnoise,'calculated faces',render='faces')
show()