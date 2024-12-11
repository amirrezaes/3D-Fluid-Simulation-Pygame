"""Utilities and physics calculations"""

from math import sqrt
from config import SimulationConfig as Config
import numpy as np
from particle import Particle

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

SURFACE_TENSION = 0.001

'''def start(
    xmin: float, xmax: float, ymin: float, zmin: float, space: float, count: int
) -> list[Particle]:
    """Creates a 3D cube of particles"""
    result = []
    x_pos, y_pos, z_pos = xmin, ymin, zmin
    particles_per_side = int(pow(count, 1/3))  # Cube root for 3D distribution
    
    for _ in range(particles_per_side):
        for _ in range(particles_per_side):
            for _ in range(particles_per_side):
                result.append(Particle(x_pos, y_pos, z_pos))
                z_pos += space
            z_pos = zmin
            y_pos += space
        y_pos = ymin
        x_pos += space
    
    return result'''


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


def calculate_density(particles: list[Particle]) -> None:
    for i, particle_1 in enumerate(particles):
        density = 0.0
        density_near = 0.0
        for particle_2 in particles[i + 1:]:
            distance = sqrt(
                (particle_1.x_pos - particle_2.x_pos) ** 2
                + (particle_1.y_pos - particle_2.y_pos) ** 2
                + (particle_1.z_pos - particle_2.z_pos) ** 2
            )
            if distance < R:
                normal_distance = 1 - distance / R
                density += normal_distance**2
                density_near += normal_distance**3
                particle_2.rho += normal_distance**2
                particle_2.rho_near += normal_distance**3
                particle_1.neighbors.append(particle_2)
        particle_1.rho += density
        particle_1.rho_near += density_near



'''def create_pressure(particles: list[Particle]) -> None:
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
        particle.z_force -= total_pressure[2]'''

def create_pressure(particles: list[Particle]) -> None:
    for particle in particles:
        press_x = 0.0
        press_y = 0.0
        press_z = 0.0
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
            normal_distance = 1 - distance / R
            total_pressure = (
                particle.press + neighbor.press
            ) * normal_distance**2 + (
                particle.press_near + neighbor.press_near
            ) * normal_distance**3
            
            pressure_vector = [
                particle_to_neighbor[0] * total_pressure / distance,
                particle_to_neighbor[1] * total_pressure / distance,
                particle_to_neighbor[2] * total_pressure / distance,
            ]
            
            neighbor.x_force += pressure_vector[0]
            neighbor.y_force += pressure_vector[1]
            neighbor.z_force += pressure_vector[2]
            press_x += pressure_vector[0]
            press_y += pressure_vector[1]
            press_z += pressure_vector[2]
            
        particle.x_force -= press_x
        particle.y_force -= press_y
        particle.z_force -= press_z


'''def calculate_viscosity(particles: list[Particle]) -> None:
    for particle in particles:
        if not particle.neighbors:
            continue
            
        # Convert positions and velocities to numpy arrays
        pos = np.array([particle.x_pos, particle.y_pos, particle.z_pos])
        vel = np.array([particle.x_vel, particle.y_vel, particle.z_vel])
        neighbor_pos = np.array([[n.x_pos, n.y_pos, n.z_pos] for n in particle.neighbors])
        neighbor_vel = np.array([[n.x_vel, n.y_vel, n.z_vel] for n in particle.neighbors])
        
        # Vectorized distance calculation
        particle_to_neighbors = neighbor_pos - pos
        distances = np.sqrt(np.sum(particle_to_neighbors**2, axis=1))
        
        # Normal vectors
        normal_vectors = particle_to_neighbors / distances[:, np.newaxis]
        relative_distances = distances / R
        
        # Velocity differences
        vel_diff = np.sum((vel - neighbor_vel) * normal_vectors, axis=1)
        
        # Mask for positive velocity differences
        positive_vel_mask = vel_diff > 0
        
        # Calculate viscosity forces
        viscosity_forces = (
            (1 - relative_distances[positive_vel_mask, np.newaxis]) * 
            SIGMA * 
            vel_diff[positive_vel_mask, np.newaxis] * 
            normal_vectors[positive_vel_mask]
        )
        
        # Update velocities
        for idx, neighbor in enumerate(particle.neighbors):
            if positive_vel_mask[idx]:
                force = viscosity_forces[np.where(positive_vel_mask)[0] == idx][0]
                # Update particle velocity
                particle.x_vel -= force[0] * 0.5
                particle.y_vel -= force[1] * 0.5
                particle.z_vel -= force[2] * 0.5
                # Update neighbor velocity
                neighbor.x_vel += force[0] * 0.5
                neighbor.y_vel += force[1] * 0.5
                neighbor.z_vel += force[2] * 0.5'''

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