from dataclasses import dataclass
from typing import Tuple
import csv
import pyvista as pv
import os

import forsys.vertex as vertex
import forsys.edge as edge
import forsys.cell as cell

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
            edges[i].gt = ground_truth_tensions[(edge_list[i][0], edge_list[i][1])]
            
        # build cell dictionary
        for i in range(len(face_list)):
            face = face_list[i]
            vertices_in_cell = []
            for v in face:
                vertices_in_cell.append(vertices[v])
            cells[i] = cell.Cell(i, vertices_in_cell)
        
    
        #for _, r in edges_temp.iterrows():
        #    edges[int(r.id)] = edge.SmallEdge(int(r.id), vertices[int(r.id1)], vertices[int(r.id2)])
        #    edges[int(r.id)].gt = round(edges_temp.loc[edges_temp['id'] == int(r.id)]['force'].iloc[0], 4)
        
        #for _, r in self.get_cells().iterrows():
        #    vlist = [edges[abs(e)].v1 if e > 0 else edges[abs(e)].v2 for e in r.edges]
        #    gt_pressure = round(r["pressures"], 4)
        #    cells[int(r.id)] = cell.Cell(int(r.id), vlist, gt_pressure=gt_pressure)
            
        
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
    
    def read_line_tensions(self):
        """
        Obtain ground truth line tensions from file.

        :return: Dictionary with vertex pairs as keys and tensions as values
        :rtype: Dict
        """
        
        line_tensions = {}
        
        # Attempt to read the line tension file
        tension_file_path = os.path.join(self.dirname, "SamplingLineTensions.csv")
        with open(tension_file_path, "r") as file:
            reader = csv.reader(file)
            edge_data = [(int(row[0]), int(row[1]), float(row[2])) for row in reader]
            
            for row in edge_data:
                line_tensions[(row[0], row[1])] = row[2]
        
        return line_tensions
        
    
    
    
    
    
