from typing import List, Tuple
import math
from dataclasses import dataclass
from particle import Particle
from numba import jit


class SpatialHashGrid:
    def __init__(self, cell_size: float):
        self.cell_size = cell_size
        self.grid: dict[Tuple[int, int, int], List[Particle]] = {}
    
    def _hash_coordinates(self, x: float, y: float, z: float) -> Tuple[int, int, int]:

        return (
            math.floor(x / self.cell_size),
            math.floor(y / self.cell_size),
            math.floor(z / self.cell_size)
        )
    
    @staticmethod
    @jit(nopython=True) 
    def _get_neighboring_cells(cell: Tuple[int, int, int]) -> List[Tuple[int, int, int]]:
        """
        Generate all potentially neighboring cell coordinates

        """
        x, y, z = cell
        neighboring_cells = [
            (x+dx, y+dy, z+dz)
            for dx in [-1, 0, 1]
            for dy in [-1, 0, 1]
            for dz in [-1, 0, 1]
        ]
        return neighboring_cells
    
    def insert(self, particle: Particle):
        """
        Insert a particle into the spatial hash grid
        """
        cell = self._hash_coordinates(particle.x_pos, particle.y_pos, particle.z_pos)
        
        if cell not in self.grid:
            self.grid[cell] = []
        
        self.grid[cell].append(particle)
    
    def clear(self):
        """
        Clear the entire spatial hash grid
        """
        self.grid.clear()
    

    def find_neighbors(
        self, 
        particle: Particle, 
        interaction_radius: float
    ) -> List[Particle]:
        """
        Find neighboring particles within interaction radius

        """
        # Reset neighbors list
        particle.neighbors = []
        
        # Get the cell of the current particle
        current_cell = self._hash_coordinates(
            particle.x_pos, 
            particle.y_pos, 
            particle.z_pos
        )
        
        # Check neighboring cells
        for cell in self._get_neighboring_cells(current_cell):
            if cell not in self.grid:
                continue
            
            # Check particles in this cell
            for potential_neighbor in self.grid[cell]:
                if potential_neighbor is particle:
                    continue
                
                # Calculate squared distance
                dx = particle.x_pos - potential_neighbor.x_pos
                dy = particle.y_pos - potential_neighbor.y_pos
                dz = particle.z_pos - potential_neighbor.z_pos
                
                squared_distance = dx*dx + dy*dy + dz*dz
                
                # Check if within interaction radius
                if squared_distance < interaction_radius * interaction_radius:
                    particle.neighbors.append(potential_neighbor)
        
        return particle.neighbors