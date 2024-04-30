from scipy.spatial import Delaunay
import numpy as np
from copy import deepcopy
from renderer import renderLineChart

class vector:
    def __init__(self,arr=[]) -> None:
        self.pos=arr[0:3] if len(arr)>=3 else None
        self.dir=arr[3:6] if len(arr)>=6 else None
        self.jacobian=np.array(arr[6:15]).reshape((3,3)) if len(arr)>=15 else None
    
    def copy(self):
        new=vector()
        new.pos=deepcopy(self.pos)
        new.dir=deepcopy(self.dir)
        new.jacobian=deepcopy(self.jacobian)
        return new

    @staticmethod
    def dist(vec1,vec2=None)->float:
        if not vec2:
            vec2=vector([0,0,0])
        return sum([(v1-v2)**2 for v1,v2 in zip(vec1.pos,vec2.pos)])**.5
    
    @staticmethod
    def difference(vec1,vec2):
        if vec1.dir and vec2.dir:
            dir=[(v1-v2) for v1,v2 in zip(vec1.dir,vec2.dir)]
            return vector([*vec1.pos,*dir])
        return vector(vec1.pos)
    
    @staticmethod
    def normed(arr):
        s=sum([(x)**2 for x in arr])**.5
        if s==0:
            return arr
        return [x/s for x in arr]
    

class Field:
    def __init__(self,file_path=None) -> None:
        self.faces=[]
        self.vertices=[]
        self.edges=[]
        if file_path is not None:
            self.parse_ply(file_path)

    def parse_ply(self,file_path):
        vertex_count = 0
        face_count = 0
        current_parser = None

        def parse_faces(line):
            parts = line.split()
            count = int(parts[0])
            if count == 3:
                indices = list(map(int, parts[1:]))
                self.faces.append(indices)

        def parse_vertices(line):
            nonlocal current_parser

            parts = line.split()
            coordinates = [float(e) for e in parts]
            self.vertices.append(vector(coordinates))
            if len(self.vertices) == vertex_count:
                current_parser=parse_faces

        def parse_header(line):
            nonlocal current_parser
            nonlocal face_count
            nonlocal vertex_count

            if line == "end_header":
                current_parser=parse_vertices
            elif line.startswith("element"):
                rest=line[8:]
                if rest[0]=='f':
                    face_count=int(rest[5:])
                elif rest[0]=='v':
                    vertex_count=int(rest[7:])

        current_parser=parse_header

        try:
            with open(file_path, 'r') as file:
                for line in file:
                    current_parser(line.strip())
        except Exception as e:
            print("error parsing ply file "+file_path+": ",e)
            data = {"error": str(e)}

        return self
    
    def calculate_faces(self):
        points=np.array([vert.pos[0:2] for vert in self.vertices])
        tri = Delaunay(points)
        self.faces=tri.simplices
        return self

    def calculate_edges(self):
        if len(self.faces)==0:
            self.calculate_faces()
        edges=[{} for _ in self.vertices]
        for v1,v2,v3 in self.faces:
            edges[v2][v1]=edges[v1][v2]=1.0/vector.dist(self.vertices[v1],self.vertices[v2])
            edges[v3][v1]=edges[v1][v3]=1.0/vector.dist(self.vertices[v1],self.vertices[v3])
            edges[v3][v2]=edges[v2][v3]=1.0/vector.dist(self.vertices[v2],self.vertices[v3])
        for v1_edges in edges:
            s=sum(v1_edges.values())
            if s!=0:
                for v2 in v1_edges: 
                    v1_edges[v2]/=s

        self.edges=edges
        return self

    def calculate_jacobian(self):
        def calculate(face,truncated=False):
            truncate=lambda x:x
            if truncated:
                truncate=lambda x:round(x,3)
            
            P=np.zeros((3,4))
            V=np.zeros((3,3))
            for i in range(3):
                P[i,3]=1
                for j in range(3):
                    P[i,j]=truncate(self.vertices[face[i]].pos[j])
            
            for i in range(3):
                for j in range(3):
                    V[i,j]=truncate(self.vertices[face[i]].dir[j])
        
            solution=np.linalg.lstsq(P,V,rcond=None)

            A_T=solution[0]
            limit=12
            if abs(A_T).max()>limit:
                if truncated: return None
                else: return calculate(face,True)

            jacobian=A_T[:3,:3].T
            return jacobian

        face_count=np.zeros((len(self.vertices)))
        def set_jacobian(v1,jacobian):
            if jacobian is None:
                return
            if self.vertices[v1].jacobian is not None:
                self.vertices[v1].jacobian+=jacobian
            else: 
                self.vertices[v1].jacobian=jacobian
            face_count[v1]+=1
        
        if len(self.faces)==0:
            self.calculate_faces()
        for face in self.faces:
            jacobian=calculate(face,True)
            for v in face:
                set_jacobian(v,jacobian)
        
        for v1 in range(len(self.vertices)):
            v=self.vertices[v1]
            if v.jacobian is not None:
                v.jacobian/=face_count[v1]
        return self

    @staticmethod
    def get_Error(calculated,expected):
        err_field=Field()
        err_field.faces=calculated.faces
        err_list=[0]*len(calculated.vertices)
        err=0
        for truth,pred in zip(calculated.vertices,expected.vertices):
            err_vector=vector.difference(truth,pred)
            cur_err=sum([x**2 for x in err_vector.dir])
            err+=cur_err
            err_list.append(cur_err)
            err_field.vertices.append(err_vector)
        err/=len(calculated.vertices)
        return err,err_field,err_list