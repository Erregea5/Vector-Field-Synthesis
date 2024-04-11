import os,sys
src_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path=os.path.dirname(src_path)+'\\data\\'
sys.path.append(src_path)

from vector_field import vector,Field
from renderer import renderField,renderLineChart,show,close
from synthesis import vectorSynthesis
import random

def main(original):
  expected=Field()
  percent_vert_present=.15
  for vertex in original.vertices:
    if random.random()<percent_vert_present:
      expected.vertices.append(vertex.copy())

  expected.calculate_jacobian()
  # renderField(expected,'original')

  def copy():
    field=Field()
    field.edges=expected.edges
    field.faces=expected.faces
    return field

  percent_present=.05
  percent_jacobian_present=.1

  sparse = copy()
  sparse_j = copy()
  sparse_jN = copy()

  for vertex in expected.vertices:
    sparse.vertices.append(vertex.copy())
    sparse_j.vertices.append(vertex.copy())
    sparse_jN.vertices.append(vertex.copy())
    sparse.vertices[-1].jacobian=None
    
    if random.random()>percent_present:
      sparse.vertices[-1].dir=None
      sparse_j.vertices[-1].dir=None
      sparse_jN.vertices[-1].dir=None
      if True:#random.random()>percent_jacobian_present:
        sparse_j.vertices[-1].jacobian=None
        sparse_jN.vertices[-1].jacobian=None

  # renderField(sparse_j,'sparse')
  # renderField(sparse_j,'sparse jacobians','jacobian')

  reconstructed=vectorSynthesis(sparse)
  reconstructed_j=vectorSynthesis(sparse_j)
  reconstructed_jN=vectorSynthesis(sparse_jN,'Nguyens')

  # renderField(reconstructed,'reconstructed without jacobian')
  # renderField(reconstructed_j,'reconstructed with jacobian')
  # renderField(reconstructed_jN,'reconstructed with Nguyen\'s method')

  _,no_err=Field.get_Error(expected,expected)
  error_field,err=Field.get_Error(reconstructed,expected)
  error_field_j,err_j=Field.get_Error(reconstructed_j,expected)
  error_field_jN,err_jN=Field.get_Error(reconstructed_jN,expected)

  print('regular error: ',err)
  print('error with jacobian: ',err_j)
  print('error with  Nguyen\'s method: ',err_jN)

  # renderField(error_field,'regular error','vertices-color')
  # renderField(error_field_j,'error with jacobian','vertices-color')
  # renderField(error_field_jN,'error with  Nguyen\'s method','vertices-color')

  # show()
  return err,err_j,err_jN


if __name__=='__main__':
  random.seed(0)
  original = Field(data_path+"bnoise.ply")
  lj,lr,ln=[],[],[]
  # main(original)
  try:
    while True: 
      r,j,n=main(original)
      lr.append(r)
      lj.append(j)
      ln.append(n)
      close()
  except KeyboardInterrupt:
    renderLineChart((lr,lj,ln),'Errors',('interp','jacobian','nguyen'))
    print('regular average:',sum(lr)/len(lr))
    print('jacobian average:',sum(lj)/len(lj))
    print('Nguyen\'s average:',sum(ln)/len(ln))
    print(len(lj),'runs ended')
  
  show()

#from results it appears that jacobian is comparable to nguyen 
#at points with values+jacobian but loses at points with just jacobian