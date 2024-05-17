import os,sys
src_path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path=os.path.dirname(src_path)+'\\data\\'
sys.path.append(src_path)
import time
import threading

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
  'divergence',
  'curl and divergence'
  ]

methods=[
  '',
  'jacobian',
  'jacobian',
  'extended jacobian',
  'extended jacobian',
  'Nguyens',
  'extended curl',
  'extended div',
  'extended curl+div'
  ]

def run_test(original,percent_present,seed,rand,debug=False):
  rand.seed(seed)
  expected=original

  def copy_field():
    field=Field()
    field.edges=expected.edges
    field.faces=expected.faces
    return field

  sparse_fields=[copy_field() for _ in range(len(titles))]
  sparse_interpolation,sparse_jacobian,sparse_jacobian_random,sparse_jacobian_augmented,sparse_jacobian_random_augmented,sparse_Nguyens,sparse_curl,sparse_divergence,sparse_curl_and_divergence = sparse_fields
  fields_with_jacobian_at_points_with_dir=[
    sparse_jacobian,
    sparse_jacobian_augmented,
    sparse_Nguyens,
    sparse_divergence,
    sparse_curl,
    sparse_curl_and_divergence
    ]
  fields_with_jacobian_randomly_placed=[
    sparse_jacobian_random,
    sparse_jacobian_random_augmented
    ]

  for vertex in expected.vertices:
    for field in sparse_fields:
      field.vertices.append(vertex.copy())
    
    if rand.random()>percent_present:
      for field in sparse_fields:
        field.vertices[-1].dir=None
      for field in fields_with_jacobian_at_points_with_dir:
        field.vertices[-1].jacobian=None

  rand.seed(seed+1)
  for i in range(len(expected.vertices)):
    if rand.random()>percent_present:
      for field in fields_with_jacobian_randomly_placed:
        field.vertices[i].jacobian=None
  
  sparse_interpolation.remove_jacobian()

  # for i in range(len(sparse_j.vertices)):
  #   r=sparse_j.vertices[i].dir is not None
  #   for j in sparse_j.edges[i]:
  #     if sparse_j.vertices[j].dir:
  #       r=True
  #   if rand.random()>percent_jacobian_present or not r:
  #     sparse_j.vertices[i].jacobian=None
  #     sparse_N.vertices[i].jacobian=None

  time_taken=[]
  reconstructed_fields:list[Field]=[]
  for field,method in zip(sparse_fields,methods):
    if method.startswith('extended '):
      method=method[9:]
      time0=time.time()
      reconstructed_fields.append(vectorSynthesisWithAugmentedMatrix(field,method))
      time_taken.append(time.time()-time0)
      continue
    time0=time.time()
    reconstructed_fields.append(vectorSynthesis(field,method,None,True))
    time_taken.append(time.time()-time0)

  curl_error_list=[]
  div_error_list=[]
  error_list=[]
  for field in reconstructed_fields:
    field.remove_jacobian()
    field.calculate_jacobian()
    curl_error_list.append(Field.get_curl_error(field,expected))
    div_error_list.append(Field.get_div_error(field,expected))
    error_list.append(Field.get_direction_error(field,expected))

  if debug:
    for message,errors,c_errors,d_errors in zip(titles,error_list,curl_error_list,div_error_list):
      print('error with '+message+': dir:',errors[0],', curl:',c_errors[0],', div:',d_errors[0])

  def render(i=''):
    renderField(expected,i+' original')
    # renderField(expected,i+' edges','faces')
    renderField(sparse_jacobian,i+' sparse')
    # renderField(sparse_jacobian,i+' sparse jacobians','jacobian')
    # for field,title in zip(reconstructed_fields,titles):
    #   renderField(field,i+' reconstructed with '+title)
    for field,err,title in zip(reconstructed_fields,error_list,titles):
      renderField(err[1],i+' error with '+title,'vertices-color')
      renderField(field,i+' reconstructed with '+title,'vertices',True)
    #   renderHistogram(err[2],i+' '+title+' error histogram')
    # for jerr,title in zip(curl_error_list,titles):
    #   renderField(jerr[1],i+' curl error with '+title,'vertices-color')
    # for jerr,title in zip(div_error_list,titles):
    #   renderField(jerr[1],i+' divergence error with '+title,'vertices-color')

    show()

  err_values=[]
  for err,cerr,derr in zip(error_list,curl_error_list,div_error_list):
    err_values.append((err[0],cerr[0],derr[0]))
  return err_values,time_taken,render

def init_field(percent_vert_present=.1):
  original = Field(data_path+"bnoise.ply")
  expected=Field()
  for vertex in original.vertices:
    if random.random()<percent_vert_present:
      expected.vertices.append(vertex.copy())

  expected.calculate_jacobian().calculate_edges()
  return expected

def main(expected):
  seed=0
  total_tests=50
  # on my machine 2 threads can run on a core so i run half of the compute on another thread
  num_extra_threads=0
  tests_per_thread=total_tests//(num_extra_threads+1)

  errors=[[0]*total_tests for _ in range(len(titles))]
  curl_errors=[[0]*total_tests for _ in range(len(titles))]
  div_errors=[[0]*total_tests for _ in range(len(titles))]
  time_taken=[[0]*total_tests for _ in range(len(titles))]

  def run_tests(start,end):
    try:
      rand=random.Random()
      for i in range(start,end): 
        percent_present=.01*(i+3)
        # print('iteration: ', i)
        
        new_errors,new_time_taken,render=run_test(expected,percent_present,seed,rand)
        for j in range(len(new_errors)):
          errors[j][i]=(new_errors[j][0])
          curl_errors[j][i]=(new_errors[j][1])
          div_errors[j][i]=(new_errors[j][2])
        for j in range(len(new_time_taken)):
          time_taken[j][i]=(new_time_taken[j])
        
    except KeyboardInterrupt:
      pass

  def show_results():
    close()
    renderLineChart(errors,'Errors',titles)
    renderLineChart(curl_errors,'curl Errors',titles)
    renderLineChart(div_errors,'divergence Errors',titles)
    renderLineChart(time_taken,'Calculation Time',titles)
    run_length=len(errors[0])
    for title,err,cerr,derr in zip(titles,errors,curl_errors,div_errors):
      print(title+' average error:',sum(err)/run_length,' curl error:',sum(cerr)/run_length,' div error:',sum(derr)/run_length)
    print(run_length,'runs ended')
    show()

  time0=time.time()
  threads=[threading.Thread(target=run_tests,args=(i*tests_per_thread,(i+1)*tests_per_thread)) for i in range(num_extra_threads)]
  for thread in threads:
    thread.start()
  run_tests((num_extra_threads)*tests_per_thread,total_tests)
  for thread in threads:
    thread.join()
  print("total time: ",(time.time()-time0))
  show_results()

if __name__ == '__main__':
  single_test=True
  random.seed(0)
  expected=init_field(.1)  
  if single_test:
    time0=time.time()
    render=run_test(expected,.2,3,random,True)[2]
    print("total time: ",(time.time()-time0))
    render()
  else:
    main(expected)
