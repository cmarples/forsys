from dataclasses import dataclass
from typing import Tuple
import csv
import pyvista as pv
import os

import forsys.vertex as vertex
import forsys.edge as edge
import forsys.cell as cell

import numpy as np


@dataclass
class Chaste:
    """Class interface with Chaste-generated files
    
    :param fname: Path to the input Chaste directory
    :type fname: str
    :param fname: Path to the Chaste file to read
    :type fname: str
    """
    dirname: str
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
        
        print("Calling Chaste reader")
        
        vertices = {}
        edges  = {}
        cells = {}
        vertex_list, edge_list, face_list = self.read_vtu_file()
        ground_truth_tensions = self.read_line_tensions()
        
        # build vertex dictionary
        for i in range(len(vertex_list)):
            vertices[i] = vertex.Vertex(i, vertex_list[i][0], vertex_list[i][1])
            
        # build edge dictionary
        for i in range(len(edge_list)):
            vertex_1 = vertices[edge_list[i][0]]
            vertex_2 = vertices[edge_list[i][1]]
            edges[i] = edge.SmallEdge(i, vertex_1, vertex_2)
            #edges[i].gt = ground_truth_tensions[(edge_list[i][0], edge_list[i][1])]
            edges[i].gt = 1.0
            
        # build cell dictionary
        for i in range(len(face_list)):
            face = face_list[i]
            vertices_in_cell = []
            for v in face:
                vertices_in_cell.append(vertices[v])
            cells[i] = cell.Cell(i, vertices_in_cell)
        
        '''
        # shuffle test (is edge ordering important)
        import random
        edge_items = list(edges.items())
        random.shuffle(edge_items)
        edges = dict(edge_items)
        '''
        
        '''
        # translation test
        xs = [v.x for v in vertices.values()]
        ys = [v.y for v in vertices.values()]
        for v in vertices.values():
            v.x = (v.x - np.mean(xs)) / (max(xs)-min(xs))
            v.y = (v.y - np.mean(ys)) / (max(ys)-min(ys))
        '''
        
        '''
        # rotation test
        for v in vertices.values():
            v.x, v.y = -v.y, v.x
        '''
        
        
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
        
        '''
        # obtain edges from vertex and face data
        edge_set = set()
        for face in face_list:
            n = len(face)
            for j in range(n):
                vertex_1 = face[j]
                vertex_2 = face[(j+1) % n] # modular arithmetic closes polygon
                edge_set.add(tuple(sorted((vertex_1, vertex_2))))

        edge_list = list(edge_set)
        '''
        
        print("edges")
        edge_set = set()
        edge_list = []
        for face in face_list:
            n = len(face)
            for j in range(n):
                e = tuple(sorted((face[j], face[(j+1)%n])))
                if e not in edge_set:
                    edge_set.add(e)
                    edge_list.append(e)
        edge_list = sorted(edge_set, key=lambda e: (min(e), max(e)))
        
        return vertex_list, edge_list, face_list
    
    def read_line_tensions(self):
        """
        Obtain ground truth line tensions from file.

        :return: Dictionary with vertex pairs as keys and tensions as values
        :rtype: Dict
        """
        
        line_tensions = {}
        
        # Attempt to read the line tension file
        tension_file_path = os.path.join(self.dirname, "GroundTruthLineTensions.csv")
        with open(tension_file_path, "r") as file:
            reader = csv.reader(file)
            edge_data = [(int(row[0]), int(row[1]), float(row[2])) for row in reader]
            
            for row in edge_data:
                line_tensions[(row[0], row[1])] = row[2]
        
        return line_tensions
        
    
    
    
    
    
