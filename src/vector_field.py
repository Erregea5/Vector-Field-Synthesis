class vector:
    def __init__(self,arr) -> None:
        self.pos=arr[0:3]
        if len(arr)==6:
            self.dir=arr[3:]
        else: self.dir=None
    
    @staticmethod
    def dist(vec1,vec2)->float:
        return sum([(v1-v2)**2 for v1,v2 in zip(vec1.pos,vec2.pos)])**.5
    

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

    def calculate_edges(self):
        edges=[{} for _ in self.vertices]
        for v1,v2,v3 in self.faces:
            edges[v2][v1]=edges[v1][v2]=vector.dist(self.vertices[v1],self.vertices[v2])
            edges[v3][v1]=edges[v1][v3]=vector.dist(self.vertices[v1],self.vertices[v3])
            edges[v3][v2]=edges[v2][v3]=vector.dist(self.vertices[v2],self.vertices[v3])
        self.edges=edges
        return self
    
