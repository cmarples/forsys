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
    
    :param dirname: Path to the input Chaste directory
    :type dirname: str
    :param fname: Path to the Chaste file to read
    :type fname: str
    :param isclip: Optional parameter, clipping the tissue boundary if true
    :type isclip: bool
    """
    dirname: str
    fname: str
    
    def __post_init__(self):
        self.vertices, self.edges, self.cells = self.create_lattice()
        self.remove_periodic_edges()
        
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
        #ground_truth_tensions = self.read_line_tensions()
        
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
    
    
    
    def remove_periodic_edges(self, periodic_in_x=True, periodic_in_y=True):
        """
        Removes all edges that wrap around a periodic boundary.
        Allows periodic boundary simulation from Chaste to be used as inputs to ForSys. 
        """
        
        # Calculate box width and height by finding maximum and minumum vertex coordinates
        xs = [v.x for v in self.vertices.values()]
        ys = [v.y for v in self.vertices.values()]
    
        box_width = max(xs) - min(xs)
        box_height = max(ys) - min(ys)
        
        
        # Detect wrapped edges
        edges_to_remove = []
        for eid, e in self.edges.items():
            dx = abs(e.v1.x - e.v2.x)
            dy = abs(e.v1.y - e.v2.y)
    
            periodic = False
            
            if periodic_in_x and dx > box_width / 2:
                periodic = True
            if periodic_in_y and dy > box_height / 2:
                periodic = True
    
            if periodic:
                edges_to_remove.append(eid)
    
        # Remove these edges
        for eid in edges_to_remove:
            del self.edges[eid]
        
    
    
    
    def clip_tissue(self):
        """
        Remove all cells that are at the boundary
        i.e. cells containing at least one vertex that has fewer than 3 neighbours.
        This is done to call ForSys on 'bulk tissue' using in-silico data with free boundaries.
        """
        
        
        edge_lookup = {}

        for eid, e in self.edges.items():
            v1_id = e.v1.id
            v2_id = e.v2.id
            key = (min(v1_id, v2_id), max(v1_id, v2_id))
            edge_lookup[key] = eid
            
        # Find the number of cells belonging to each edge
        edge_cell_count = {eid: 0 for eid in self.edges}
        for cid, c in self.cells.items():
            verts = c.vertices
            n = len(verts)
    
            for i in range(n):
                v1_id = verts[i].id
                v2_id = verts[(i + 1) % n].id
                key = (min(v1_id, v2_id), max(v1_id, v2_id))
    
                eid = edge_lookup[key]
                edge_cell_count[eid] += 1
                
        
        # Find the boundary vertices
        boundary_vertex_ids = set()
        for eid, count in edge_cell_count.items():
            if count == 1:  # boundary edge
                e = self.edges[eid]
                boundary_vertex_ids.add(e.v1.id)
                boundary_vertex_ids.add(e.v2.id)
            
        # Keep interior cells
        kept_cells = {}
        for cid, c in self.cells.items():
            if not any(v.id in boundary_vertex_ids for v in c.vertices):
                kept_cells[cid] = c   
        
        # Keep interior edges
        used_edge_ids = set()
        for cid, c in kept_cells.items():
            verts = c.vertices
            n = len(verts)
    
            for i in range(n):
                v1_id = verts[i].id
                v2_id = verts[(i + 1) % n].id
                key = (min(v1_id, v2_id), max(v1_id, v2_id))
    
                eid = edge_lookup[key]
                used_edge_ids.add(eid)
    
        kept_edges = {eid: self.edges[eid] for eid in used_edge_ids}    
        
        # Keep interior vertices
        used_vertex_ids = set()
        for eid in used_edge_ids:
            e = self.edges[eid]
            used_vertex_ids.add(e.v1.id)
            used_vertex_ids.add(e.v2.id)
    
        kept_vertices = {vid: self.vertices[vid] for vid in used_vertex_ids}
            
            
        self.vertices = kept_vertices
        self.edges = kept_edges
        self.cells = kept_cells
    
    
    
    
