"""Utilities and physics calculations"""

from math import sqrt
from config import SimulationConfig as Config
import numpy as np
from particle import Particle
from hashg import SpatialHashGrid
from typing import List

(
    N,
    SIM_W,
    BOTTOM,
    SIM_DEPTH,
    DAM,
    DAM_BREAK,
    G,
    SPACING,
    K,
    K_NEAR,
    REST_DENSITY,
    R,
    SIGMA,
    MAX_VEL,
    WALL_DAMP,
    VEL_DAMP,
) = Config().get_simulation_parameters()


def start(xmin: float, xmax: float, ymin: float, zmin: float, space: float, count: int) -> list[Particle]:
    """Creates a 3D cube of particles using vectorized operations"""
    particles_per_side = int(pow(count, 1/3))  # Cube root for 3D distribution
    
    # Create coordinate grids
    x = np.linspace(xmin, xmin + (particles_per_side-1)*space, particles_per_side)
    y = np.linspace(ymin, ymin + (particles_per_side-1)*space, particles_per_side)
    z = np.linspace(zmin, zmin + (particles_per_side-1)*space, particles_per_side)
    
    # Generate all combinations of coordinates
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # Flatten coordinates and create particles
    positions = np.stack([X.flatten(), Y.flatten(), Z.flatten()], axis=1)
    
    return [Particle(pos[0], pos[1], pos[2]) for pos in positions]


def calculate_density( particles: list[Particle], interaction_radius: float = R) -> None:
    grid = SpatialHashGrid(cell_size=interaction_radius)
    
    # Clear and insert particles
    grid.clear()
    for particle in particles:
        # Reset particle properties
        particle.rho = 0.0
        particle.rho_near = 0.0
        particle.neighbors = []
        
        # Insert into grid
        grid.insert(particle)
    
    # Calculate densities
    for particle in particles:
        # Find neighbors efficiently
        neighbors = grid.find_neighbors(particle, interaction_radius)
        
        # Density calculation (similar to original method)
        for neighbor in neighbors:
            distance = sqrt(
                (particle.x_pos - neighbor.x_pos)**2 +
                (particle.y_pos - neighbor.y_pos)**2 +
                (particle.z_pos - neighbor.z_pos)**2
            )
            
            if distance < interaction_radius:
                normal_distance = 1 - distance / interaction_radius
                
                density_contribution = normal_distance**2
                density_near_contribution = normal_distance**3
                
                particle.rho += density_contribution
                particle.rho_near += density_near_contribution



def create_pressure(particles: list[Particle]) -> None:
    for particle in particles:
        if not particle.neighbors:
            continue
            
        # Convert positions to numpy arrays
        pos = np.array([particle.x_pos, particle.y_pos, particle.z_pos])
        neighbor_pos = np.array([[n.x_pos, n.y_pos, n.z_pos] for n in particle.neighbors])
        
        # Vectorized distance calculation
        particle_to_neighbors = neighbor_pos - pos
        distances = np.sqrt(np.sum(particle_to_neighbors**2, axis=1))
        
        # Pressure calculations
        normal_distances = 1 - distances / R
        neighbor_press = np.array([n.press for n in particle.neighbors])
        neighbor_press_near = np.array([n.press_near for n in particle.neighbors])
        
        total_pressures = (
            (particle.press + neighbor_press) * normal_distances**2 +
            (particle.press_near + neighbor_press_near) * normal_distances**3
        )
        
        # Pressure vectors calculation
        pressure_vectors = (particle_to_neighbors.T * (total_pressures / distances)).T
        
        # Update forces
        for idx, neighbor in enumerate(particle.neighbors):
            neighbor.x_force += pressure_vectors[idx, 0]
            neighbor.y_force += pressure_vectors[idx, 1]
            neighbor.z_force += pressure_vectors[idx, 2]
        
        # Update particle forces
        total_pressure = np.sum(pressure_vectors, axis=0)
        particle.x_force -= total_pressure[0]
        particle.y_force -= total_pressure[1]
        particle.z_force -= total_pressure[2]


def calculate_viscosity(particles: list[Particle]) -> None:
    for particle in particles:
        for neighbor in particle.neighbors:
            particle_to_neighbor = [
                neighbor.x_pos - particle.x_pos,
                neighbor.y_pos - particle.y_pos,
                neighbor.z_pos - particle.z_pos,
            ]
            distance = sqrt(
                particle_to_neighbor[0] ** 2 
                + particle_to_neighbor[1] ** 2 
                + particle_to_neighbor[2] ** 2
            )
            normal_p_to_n = [
                particle_to_neighbor[0] / distance,
                particle_to_neighbor[1] / distance,
                particle_to_neighbor[2] / distance,
            ]
            relative_distance = distance / R
            velocity_difference = (
                (particle.x_vel - neighbor.x_vel) * normal_p_to_n[0]
                + (particle.y_vel - neighbor.y_vel) * normal_p_to_n[1]
                + (particle.z_vel - neighbor.z_vel) * normal_p_to_n[2]
            )
            
            if velocity_difference > 0:
                viscosity_force = [
                    (1 - relative_distance) * SIGMA * velocity_difference * normal_p_to_n[0],
                    (1 - relative_distance) * SIGMA * velocity_difference * normal_p_to_n[1],
                    (1 - relative_distance) * SIGMA * velocity_difference * normal_p_to_n[2],
                ]
                particle.x_vel -= viscosity_force[0] * 0.5
                particle.y_vel -= viscosity_force[1] * 0.5
                particle.z_vel -= viscosity_force[2] * 0.5
                neighbor.x_vel += viscosity_force[0] * 0.5
                neighbor.y_vel += viscosity_force[1] * 0.5
                neighbor.z_vel += viscosity_force[2] * 0.5