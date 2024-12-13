"""Configuration module for fluid dynamics particle simulation parameters."""

# Simulation domain parameters
NUM_PARTICLES = 300       # not exact, but close
SIMULATION_WIDTH = 0.35  # Horizontal extent of the simulation space (x-axis)
SIMULATION_DEPTH = 0.12  # Vertical extent of the simulation space (z-axis)
GROUND_LEVEL = 0  # Base y-coordinate of the simulation space
DAM_POSITION = -0.13  # Spatial location of the dam within the simulation
DAM_BREAK_FRAME = 200  # Simulation frame at which dam rupture occurs

# Physics simulation parameters
GRAVITY_ACCELERATION = 0.02 * 0.25  # Gravitational acceleration applied to particles
PARTICLE_SPACING = 0.085  # Interparticle distance used for pressure calculations
PRESSURE_COEFFICIENT = PARTICLE_SPACING / 1000.0  # Baseline pressure scaling factor
NEAR_PRESSURE_COEFFICIENT = PRESSURE_COEFFICIENT * 100  # Pressure scaling for closely positioned particles
REFERENCE_PARTICLE_DENSITY = 3  # Baseline particle density for comparative calculations
INTERACTION_RADIUS = PARTICLE_SPACING * 1.2  # Maximum distance for particle neighborhood interactions
VISCOSITY_DAMPING = 0.3  # Coefficient controlling fluid viscosity and energy dissipation
MAXIMUM_PARTICLE_VELOCITY = 1.0  # Velocity threshold to prevent simulation instability
WALL_INTERACTION_DAMPING = 0.07  # Coefficient for constraining particles within simulation boundaries
VELOCITY_REDUCTION_FACTOR = 0.5  # Scaling factor for velocity mitigation


class SimulationConfig:
    """Manages configuration parameters for the particle-based fluid dynamics simulation."""

    def __init__(self):
        """Initialize the configuration parameters."""
        return None

    def get_simulation_parameters(self):
        """
        Retrieves comprehensive simulation and physics configuration parameters.
        
        Returns:
            tuple: Ordered configuration parameters for simulation setup
        """
        return (
            NUM_PARTICLES,
            SIMULATION_WIDTH,
            GROUND_LEVEL,
            SIMULATION_DEPTH,
            DAM_POSITION,
            DAM_BREAK_FRAME,
            GRAVITY_ACCELERATION,
            PARTICLE_SPACING,
            PRESSURE_COEFFICIENT,
            NEAR_PRESSURE_COEFFICIENT,
            REFERENCE_PARTICLE_DENSITY,
            INTERACTION_RADIUS,
            VISCOSITY_DAMPING,
            MAXIMUM_PARTICLE_VELOCITY,
            WALL_INTERACTION_DAMPING,
            VELOCITY_REDUCTION_FACTOR
        )