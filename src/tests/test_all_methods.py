import os,sys
src_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path=os.path.dirname(src_path)+'\\data\\'
sys.path.append(src_path)
from matplotlib import pyplot as plt

from vector_field import vector,Field
from renderer import renderField,renderLineChart,renderHistogram,show,close
from synthesis import vectorSynthesis
import random

def main(original,e):
  expected=Field()
  percent_vert_present=.15
  for vertex in original.vertices:
    if random.random()<percent_vert_present:
      expected.vertices.append(vertex.copy())

  expected.calculate_jacobian()
  expected.calculate_edges()

  def copy_field():
    field=Field()
    field.edges=expected.edges
    field.faces=expected.faces
    return field

  percent_present=e
  percent_jacobian_present=e#.05

  sparse = copy_field()
  sparse_j = copy_field()
  sparse_jN = copy_field()

  for vertex in expected.vertices:
    sparse.vertices.append(vertex.copy())
    sparse_j.vertices.append(vertex.copy())
    sparse_jN.vertices.append(vertex.copy())
    sparse.vertices[-1].jacobian=None
    
    if random.random()>percent_present:
      sparse.vertices[-1].dir=None
      sparse_j.vertices[-1].dir=None
      sparse_jN.vertices[-1].dir=None
    
      #if True only points with direction have jacobian, 
      #otherwise some points without direction also have jacobian
      #if indented and True no points have jacobian
      #if not true and indented its random
    # if random.random()>percent_jacobian_present:
      sparse_j.vertices[-1].jacobian=None
      sparse_jN.vertices[-1].jacobian=None

  # for i in range(len(sparse_j.vertices)):
  #   dist=100
  #   r=sparse_j.vertices[i].dir is not None
  #   # for j in range(len(sparse_j.vertices)):
  #   #   if sparse_j.vertices[j].dir:
  #   #     r=True
  #   #     dist_=sum([(x-y)**2 for (x,y) in zip(sparse_j.vertices[i].pos,sparse_j.vertices[j].pos)])**.5
  #   #     if dist_<dist:
  #   #       dist=dist_
  #   # if not r:
  #   #   print("HELPPPPP")
  #   for j in sparse_j.edges[i]:
  #     if sparse_j.vertices[j].dir:
  #       r=True
  #   if random.random()>percent_jacobian_present or not r:
  #     sparse_j.vertices[i].jacobian=None
  #     sparse_jN.vertices[i].jacobian=None

  reconstructed=vectorSynthesis(sparse)
  reconstructed_j=vectorSynthesis(sparse_j,'merged',expected)
  reconstructed_jN=vectorSynthesis(sparse_jN,'Nguyens')

  error_field,err,err_list=Field.get_Error(reconstructed,expected)
  error_field_j,err_j,err_list_j=Field.get_Error(reconstructed_j,expected)
  error_field_jN,err_jN,err_list_jN=Field.get_Error(reconstructed_jN,expected)

  print('regular error: ',err)
  print('error with jacobian: ',err_j)
  print('error with  Nguyen\'s method: ',err_jN)

  def render(i=''):
    renderField(expected,i+' original')
    renderField(expected,i+' edges','faces')
    renderField(sparse_j,i+' sparse')
    renderField(sparse_j,i+' sparse jacobians','jacobian')
    renderField(reconstructed,i+' reconstructed without jacobian')
    renderField(reconstructed_j,i+' reconstructed with jacobian')
    renderField(reconstructed_jN,i+' reconstructed with Nguyen\'s method')
    # renderField(error_field,i+' regular error','vertices-color')
    renderField(error_field_j,i+' error with jacobian','vertices-color')
    renderField(error_field_jN,i+' error with  Nguyen\'s method','vertices-color')
    renderHistogram(err_list_j,i+' jacobian error histogram')
    renderHistogram(err_list_jN,i+' Nguyen\'s error histogram')
    show()
  # render()
  return err,err_j,err_jN,render


if __name__=='__main__':
  random.seed(0)
  original = Field(data_path+"bnoise.ply")
  lj,lr,ln=[],[],[]
  # random.seed(1)
  # main(original,.15)[3]()
  # exit()
  i=0
  try:
    while True: 
      i+=1
      random.seed(1)
      percent_present=.02+.03*i
      print('iteration: ', i)
      r,j,n,render=main(original,percent_present)
      lr.append(r)
      lj.append(j)
      ln.append(n)
      # if r<j or r<n:
      #   render(str(i))
      
  except KeyboardInterrupt:
    close()
    renderLineChart((lr,lj,ln),'Errors',('interp','jacobian','nguyen'))
    print('regular average:',sum(lr)/len(lr))
    print('jacobian average:',sum(lj)/len(lj))
    print('Nguyen\'s average:',sum(ln)/len(ln))
    print(len(lj),'runs ended')
  show()

# show methods accross different 

#from results it appears that jacobian is comparable to nguyen 
#at points with values+jacobian but loses at points with just jacobian

#debug singular matrix
#find why interpolation gets the same results sometimes
#more buckets in color plots
#indexed seeding
#histogram of erro
#n=m info
#lic


"""
different controls:
  I can use only points with or neighboring points with direction | random
  I can use one | all neighbors in jacobian calculation
  I can always | sometimes truncate

  -too much truncation reduces usefulness of jacobian
  -what are the chances user only chooses  
  -when having values 
"""