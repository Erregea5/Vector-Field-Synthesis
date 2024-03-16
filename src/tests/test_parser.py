import os,sys
src_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path=os.path.dirname(src_path)+'\\data\\'
sys.path.append(src_path)

from vector_field import Field
from renderer import renderField,show

bnoise = Field(data_path+"bnoise.ply")

renderField(bnoise,'0')
renderField(bnoise,'1')
show()