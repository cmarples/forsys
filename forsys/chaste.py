from dataclasses import dataclass
import pyvista as pv
from typing import Tuple

import forsys.vertex as vertex
import forsys.edge as edge
import forsys.cell as cell

@dataclass
class Chaste:
    """Class interface with Chaste-generated files
    
    :param fname: Path to the Chaste file to read
    :type fname: str
    """
    fname: str
    
    def __post_init__(self):
        self.vertices, self.edges, self.cells = self.create_lattice()
        
    def create_lattice(self) -> Tuple:
        """
        Create vertices, edges and cells from an output Chaste .vtu file. 
        All necessary steps are taken by this call.

        :return: Three dictionaries with the vertices, edges and cells respectively   
        :rtype: Tuple
        """
        
        vertices = {}
        cells = {}
        edges  = {}
        vertex_list, edge_list, face_list = self.read_vtu_file()
        
        # build vertex dictionary
        for i in range(len(vertex_list)):
            vertices[i] = vertex.Vertex(i, vertex_list[i][0], vertex_list[i][1])
        
            
        
        return vertices, edges, cells
    
    def read_vtu_file(self):
        # Read file
        mesh = pv.read(self.fname)
        vertex_list = mesh.points.tolist()
        face_array = mesh.cells # this is a 1D array
        
        # convert face_array to a list of individual faces
        face_list = []
        index = 0
        while index < len(face_array):
            n = face_array[index]
            vertices_in_face = face_array[index+1 : index+1+n]
            face_list.append(vertices_in_face)
            index += n + 1
        
        # obtain edges from vertex and face data
        edge_set = set()
        for face in face_list:
            n = len(face)
            for j in range(n):
                vertex_1 = face[j]
                vertex_2 = face[(j+1) % n] # modular arithmetic closes polygon
                edge_set.add(tuple(sorted((vertex_1, vertex_2))))

        edge_list = list(edge_set)
        
        return vertex_list, edge_list, face_list
    

        
    
    
    
    
    
