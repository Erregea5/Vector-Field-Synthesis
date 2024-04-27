import os,sys
src_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path=os.path.dirname(src_path)+'\\data\\'
sys.path.append(src_path)
import time

from vector_field import vector,Field
from renderer import renderField,renderLineChart,renderHistogram,show,close
from synthesis import vectorSynthesis
import random

def main(original,percent_present,seed):
  expected=original

  def copy_field():
    field=Field()
    field.edges=expected.edges
    field.faces=expected.faces
    return field

  sparse = copy_field()
  sparse_j = copy_field()
  sparse_jr = copy_field()
  sparse_N = copy_field()

  for vertex in expected.vertices:
    sparse.vertices.append(vertex.copy())
    sparse_j.vertices.append(vertex.copy())
    sparse_jr.vertices.append(vertex.copy())
    sparse_N.vertices.append(vertex.copy())
    sparse.vertices[-1].jacobian=None
    
    if random.random()>percent_present:
      sparse.vertices[-1].dir=None
      sparse_j.vertices[-1].dir=None
      sparse_jr.vertices[-1].dir=None
      sparse_N.vertices[-1].dir=None

      sparse_j.vertices[-1].jacobian=None
      sparse_N.vertices[-1].jacobian=None

  random.seed(seed+1)
  for i in range(len(expected.vertices)):
    if random.random()>percent_present:
      sparse_jr.vertices[i].jacobian=None

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
  #     sparse_N.vertices[i].jacobian=None

  time0=time.time()
  reconstructed=vectorSynthesis(sparse)
  dtime=time.time()-time0
  time0=time.time()
  reconstructed_j=vectorSynthesis(sparse_j,'merged')
  dtime_j=time.time()-time0
  time0=time.time()
  reconstructed_jr=vectorSynthesis(sparse_jr,'merged')
  dtime_jr=time.time()-time0
  time0=time.time()
  reconstructed_N=vectorSynthesis(sparse_N,'Nguyens')
  dtime_N=time.time()-time0

  error_field,err,err_list=Field.get_Error(reconstructed,expected)
  error_field_j,err_j,err_list_j=Field.get_Error(reconstructed_j,expected)
  error_field_jr,err_jr,err_list_jr=Field.get_Error(reconstructed_jr,expected)
  error_field_N,err_N,err_list_N=Field.get_Error(reconstructed_N,expected)

  print('regular error: ',err)
  print('error with jacobian: ',err_j)
  print('error with jacobian at random points: ',err_jr)
  print('error with  Nguyen\'s method: ',err_N)

  def render(i=''):
    renderField(expected,i+' original')
    renderField(expected,i+' edges','faces')
    renderField(sparse_j,i+' sparse')
    renderField(sparse_j,i+' sparse jacobians','jacobian')
    renderField(reconstructed,i+' reconstructed without jacobian')
    renderField(reconstructed_j,i+' reconstructed with jacobian')
    renderField(reconstructed_N,i+' reconstructed with Nguyen\'s method')
    # renderField(error_field,i+' regular error','vertices-color')
    renderField(error_field_j,i+' error with jacobian','vertices-color')
    renderField(error_field_N,i+' error with  Nguyen\'s method','vertices-color')
    renderHistogram(err_list_j,i+' jacobian error histogram')
    renderHistogram(err_list_N,i+' Nguyen\'s error histogram')
    show()
  # render()
  return [err,err_j,err_jr,err_N,dtime,dtime_j,dtime_jr,dtime_N],render


if __name__=='__main__':
  random.seed(0)
  original = Field(data_path+"bnoise.ply")
  expected=Field()
  percent_vert_present=.15
  for vertex in original.vertices:
    if random.random()<percent_vert_present:
      expected.vertices.append(vertex.copy())

  expected.calculate_jacobian()
  expected.calculate_edges()

  data=[[] for _ in range(8)]
  # random.seed(1)
  # main(original,.15)[6]()
  # exit()
  i=0
  try:
    while True: 
      i+=1
      random.seed(3)
      percent_present=.01*i
      print('iteration: ', i)
      newdata,render=main(expected,percent_present,3)
      for j in range(len(newdata)):
        data[j].append(newdata[j])
      
  except KeyboardInterrupt:
    close()
    labels=('interp','jacobian','jacobian at random points','nguyen')
    renderLineChart(data[:4],'Errors',labels)
    renderLineChart(data[4:],'Calculation Time',labels)
    run_length=len(data[0])
    print('regular average:',sum(data[0])/run_length)
    print('jacobian average:',sum(data[1])/run_length)
    print('jacobian at random points average:',sum(data[2])/run_length)
    print('Nguyen\'s average:',sum(data[3])/run_length)
    print(run_length,'runs ended')
  show()

# try on new data sets and smaller intervals
# fix nguyen's method
# write methodology of each system
# add curl and divergence

"""
different controls:
  I can use only points with or neighboring points with direction | random
  I can use one | all neighbors in jacobian calculation
  I can always | sometimes truncate

  -too much truncation reduces usefulness of jacobian
  -what are the chances user only chooses  
  -when having values 
"""