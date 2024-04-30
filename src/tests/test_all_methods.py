import os,sys
src_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path=os.path.dirname(src_path)+'\\data\\'
sys.path.append(src_path)
import time

from vector_field import vector,Field
from renderer import renderField,renderLineChart,renderHistogram,show,close
from synthesis import vectorSynthesis,vectorSynthesisWithAugmentedMatrix
import random


titles=[
  'interpolation',
  'jacobian at points with direction',
  'jacobian at random points',
  'jacobian at points with direction with extended matrix',
  'jacobian at random points with extended matrix',
  'Nguyen\'s method',
  'curl',
  'divergence'
  ]

def main(original,percent_present,seed):
  random.seed(seed)
  expected=original

  def copy_field():
    field=Field()
    field.edges=expected.edges
    field.faces=expected.faces
    return field

  sparse_fields=[copy_field() for _ in range(8)]
  sparse_interpolation,sparse_jacobian,sparse_jacobian_random,sparse_jacobian_augmented,sparse_jacobian_random_augmented,sparse_Nguyens,sparse_curl,sparse_divergence = sparse_fields
  fields_with_jacobian_at_points_with_dir=[
    sparse_jacobian,
    sparse_jacobian_augmented,
    sparse_Nguyens,
    sparse_divergence,
    sparse_curl
    ]
  fields_with_jacobian_randomly_placed=[
    sparse_jacobian_random,
    sparse_jacobian_random_augmented
    ]

  for vertex in expected.vertices:
    for field in sparse_fields:
      field.vertices.append(vertex.copy())
    sparse_interpolation.vertices[-1].jacobian=None
    
    if random.random()>percent_present:
      for field in sparse_fields:
        field.vertices[-1].dir=None
      for field in fields_with_jacobian_at_points_with_dir:
        field.vertices[-1].jacobian=None

  random.seed(seed+1)
  for i in range(len(expected.vertices)):
    if random.random()>percent_present:
      for field in fields_with_jacobian_randomly_placed:
        field.vertices[i].jacobian=None

  # for i in range(len(sparse_j.vertices)):
  #   r=sparse_j.vertices[i].dir is not None
  #   for j in sparse_j.edges[i]:
  #     if sparse_j.vertices[j].dir:
  #       r=True
  #   if random.random()>percent_jacobian_present or not r:
  #     sparse_j.vertices[i].jacobian=None
  #     sparse_N.vertices[i].jacobian=None

  time_taken=[]
  reconstructed_fields=[]
  methods=['','merged','merged','jacobian','jacobian','Nguyens','curl','div']
  for field,method in zip(sparse_fields,methods):
    if method=='curl' or method=='div' or method=='jacobian':
      time0=time.time()
      reconstructed_fields.append(vectorSynthesisWithAugmentedMatrix(field,method))
      time_taken.append(time.time()-time0)
      continue
    time0=time.time()
    reconstructed_fields.append(vectorSynthesis(field,method))
    time_taken.append(time.time()-time0)

  error_list=[]
  for field in reconstructed_fields:
    error_list.append(Field.get_Error(field,expected))

  for message,errors in zip(titles,error_list):
    print('error with '+message+': ',errors[0])


  def render(i=''):
    renderField(expected,i+' original')
    renderField(expected,i+' edges','faces')
    renderField(sparse_jacobian,i+' sparse')
    renderField(sparse_jacobian,i+' sparse jacobians','jacobian')
    for field,title in zip(reconstructed_fields,titles):
      renderField(field,i+' reconstructed with '+title)
    # for err,title in zip(error_list,titles):
    #   renderField(err[1],i+' error with '+title,'vertices-color')
    # for err,title in zip(error_list,titles):
    #   renderHistogram(err[2],i+' '+title+' error histogram')
    show()

  err_values=[]
  for err in error_list:
    err_values.append(err[0])
  return err_values,time_taken,render


if __name__=='__main__':
  random.seed(0)
  original = Field(data_path+"bnoise.ply")
  expected=Field()
  percent_vert_present=.1
  for vertex in original.vertices:
    if random.random()<percent_vert_present:
      expected.vertices.append(vertex.copy())

  expected.calculate_jacobian()
  expected.calculate_edges()

  # main(expected,1,3)[2]()
  # exit()

  i=0
  errors=[[] for _ in range(8)]
  time_taken=[[] for _ in range(8)]
  try:
    while True: 
      i+=1
      percent_present=.01*i
      print('iteration: ', i)
      new_errors,new_time_taken,render=main(expected,percent_present,3)
      for j in range(len(new_errors)):
        errors[j].append(new_errors[j])
      for j in range(len(new_time_taken)):
        time_taken[j].append(new_time_taken[j])
      
  except KeyboardInterrupt:
    close()
    renderLineChart(errors,'Errors',titles)
    renderLineChart(time_taken,'Calculation Time',titles)
    run_length=len(errors[0])
    for title,err in zip(titles,errors):
      print(title+' average:',sum(err)/run_length)
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