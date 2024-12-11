"""Defines the Particle class."""

from math import sqrt
from config import SimulationConfig as Config

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

class Particle:
    def __init__(self, x_pos: float, y_pos: float, z_pos: float):
        self.x_pos = x_pos
        self.y_pos = y_pos
        self.z_pos = z_pos
        self.previous_x_pos = x_pos
        self.previous_y_pos = y_pos
        self.previous_z_pos = z_pos
        self.visual_x_pos = x_pos
        self.visual_y_pos = y_pos
        self.visual_z_pos = z_pos
        self.rho = 0.0
        self.rho_near = 0.0
        self.press = 0.0
        self.press_near = 0.0
        self.neighbors = []
        self.x_vel = 0.0
        self.y_vel = 0.0
        self.z_vel = 0.0
        self.x_force = 0.0
        self.y_force = -G
        self.z_force = 0.0


    def update_state(self, dam: bool):
        """
        Updates the state of the particle using leapfrog integration
        """
        dt = 1.0  # timestep
        half_dt = dt / 2.0

        # Store current position
        self.previous_x_pos = self.x_pos
        self.previous_y_pos = self.y_pos

        # First half of velocity update (v(t + dt/2) = v(t) + a(t) * dt/2)
        self.x_vel += self.x_force * half_dt
        self.y_vel += self.y_force * half_dt

        # Position update using half-step velocity (x(t + dt) = x(t) + v(t + dt/2) * dt)
        self.x_pos += self.x_vel * dt
        self.y_pos += self.y_vel * dt

        # Set visual position
        self.visual_x_pos = self.x_pos
        self.visual_y_pos = self.y_pos

        # Reset forces except gravity
        self.x_force = 0.0
        self.y_force = -G

        # Apply wall constraints and update forces
        if self.x_pos < -SIM_W:
            self.x_force -= (self.x_pos - -SIM_W) * WALL_DAMP
            self.visual_x_pos = -SIM_W

        if dam is True and self.x_pos > DAM:
            self.x_force -= (self.x_pos - DAM) * WALL_DAMP**(abs(self.x_pos - DAM))

        if self.x_pos > SIM_W:
            self.x_force -= (self.x_pos - SIM_W) * WALL_DAMP
            self.visual_x_pos = SIM_W

        if self.y_pos < BOTTOM:
            self.y_force -= (self.y_pos - SIM_W) * WALL_DAMP
            self.visual_y_pos = BOTTOM
        
        # Second half of velocity update (v(t + dt) = v(t + dt/2) + a(t + dt) * dt/2)
        self.x_vel += self.x_force * half_dt
        self.y_vel += self.y_force * half_dt

        # Calculate velocity magnitude
        velocity = sqrt(self.x_vel**2 + self.y_vel**2)

        # Apply velocity damping if needed
        if velocity > MAX_VEL:
            self.x_vel *= VEL_DAMP
            self.y_vel *= VEL_DAMP

        # Reset density and neighbors for next iteration
        self.rho = 0.0
        self.rho_near = 0.0
        self.neighbors = []

    def calculate_pressure(self):
        self.press = K * (self.rho - REST_DENSITY)
        self.press_near = K_NEAR * self.rho_near