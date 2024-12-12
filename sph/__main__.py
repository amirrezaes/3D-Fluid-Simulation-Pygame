import pygame
import numpy as np
from OpenGL.GL import *
from OpenGL.GLU import *
import psutil

from math import sqrt
from config import SimulationConfig as Config
from particle import Particle
from physics import start, calculate_density, create_pressure, calculate_viscosity

WINDOW_SIZE = (800, 600)
PARTICLE_SIZE = 50
FOV = 45.0
NEAR = 0.1
FAR = 100.0
ASPECT_RATIO = WINDOW_SIZE[0] / WINDOW_SIZE[1]

class AnimatedScatter():
    def __init__(self):
        super(AnimatedScatter, self).__init__()
        (
            self.NUM_PARTICLES, self.SIMULATION_WIDTH, self.GROUND_LEVEL, self.SIMULATION_DEPTH, 
            self.DAM_POSITION, self.DAM_BREAK_FRAME, self.GRAVITY_ACCELERATION,
            self.PARTICLE_SPACING, self.PRESSURE_COEFFICIENT, self.NEAR_PRESSURE_COEFFICIENT, 
            self.REFERENCE_PARTICLE_DENSITY, self.INTERACTION_RADIUS, self.VISCOSITY_DAMPING,
            self.MAXIMUM_PARTICLE_VELOCITY, self.WALL_INTERACTION_DAMPING, self.VELOCITY_REDUCTION_FACTOR
        ) = Config().get_simulation_parameters()

        self.simulation_state = start(
            -self.SIMULATION_WIDTH, self.DAM_POSITION, self.GROUND_LEVEL, -self.SIMULATION_DEPTH, 
            0.03, self.NUM_PARTICLES
        )
        
        self.dam_built = True
        self.frame = 0
        
        self.camera_distance = 2
        self.camera_x = 0.0
        self.camera_y = 0.0
        self.camera_z = self.camera_distance
        self.rot_x = 25
        self.rot_y = -5

        self.color_func = lambda x: glColor4f(0.0, 0.5, 0.8, 0.8)

        
        pygame.init()
        self.font = pygame.font.SysFont('Arial', 18)
        self.screen = pygame.display.set_mode(WINDOW_SIZE, pygame.DOUBLEBUF | pygame.OPENGL)
        pygame.display.set_caption('3D SPH Fluid Simulation')
        self.clock = pygame.time.Clock()
        self.running = True
        self.setup_opengl()

    def setup_opengl(self):
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_COLOR_MATERIAL)
        glEnable(GL_POINT_SMOOTH)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        
        glLightfv(GL_LIGHT0, GL_POSITION, (0, 1, 1, 0))
        glLightfv(GL_LIGHT0, GL_AMBIENT, (0.2, 0.2, 0.2, 1.0))
        glLightfv(GL_LIGHT0, GL_DIFFUSE, (0.5, 0.5, 0.5, 1.0))
        
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(FOV, ASPECT_RATIO, NEAR, FAR)
        
        glMatrixMode(GL_MODELVIEW)
        glPointSize(PARTICLE_SIZE)

    def handle_input(self):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                self.running = False
                return
            
        keys = pygame.key.get_pressed()
        
        if keys[pygame.K_LEFT]:
            self.rot_y -= 1
        if keys[pygame.K_RIGHT]:
            self.rot_y += 1
        if keys[pygame.K_UP]:
            self.rot_x -= 1
        if keys[pygame.K_DOWN]:
            self.rot_x += 1
            
        if keys[pygame.K_q]:
            self.camera_distance -= 0.1
        if keys[pygame.K_e]:
            self.camera_distance += 0.1
        
        if keys[pygame.K_v]:
            self.color_func = self.color_velocity
        
        if keys[pygame.K_p]:
            self.color_func = self.color_press
        
        if keys[pygame.K_ESCAPE]:
            self.running = False
            
    
    def draw_3d_walls(self):
        glBegin(GL_QUADS)
        glColor3f(0.0, 1.0, 0.0)

        # bottom
        glVertex3f(-self.SIMULATION_WIDTH - 0.02, self.GROUND_LEVEL, -self.SIMULATION_DEPTH)
        glVertex3f(self.SIMULATION_WIDTH  + 0.02, self.GROUND_LEVEL, -self.SIMULATION_DEPTH)
        glVertex3f(self.SIMULATION_WIDTH  + 0.02, self.GROUND_LEVEL, self.SIMULATION_DEPTH)
        glVertex3f(-self.SIMULATION_WIDTH - 0.02, self.GROUND_LEVEL, self.SIMULATION_DEPTH)

        # left
        glVertex3f(-self.SIMULATION_WIDTH - 0.02, self.GROUND_LEVEL, -self.SIMULATION_DEPTH)
        glVertex3f(-self.SIMULATION_WIDTH - 0.02, self.GROUND_LEVEL, self.SIMULATION_DEPTH)
        glVertex3f(-self.SIMULATION_WIDTH - 0.02, self.SIMULATION_WIDTH, self.SIMULATION_DEPTH)
        glVertex3f(-self.SIMULATION_WIDTH - 0.02, self.SIMULATION_WIDTH, -self.SIMULATION_DEPTH)

        # right
        glVertex3f(self.SIMULATION_WIDTH + 0.02, self.GROUND_LEVEL, -self.SIMULATION_DEPTH)
        glVertex3f(self.SIMULATION_WIDTH + 0.02, self.GROUND_LEVEL, self.SIMULATION_DEPTH)
        glVertex3f(self.SIMULATION_WIDTH + 0.02, self.SIMULATION_WIDTH, self.SIMULATION_DEPTH)
        glVertex3f(self.SIMULATION_WIDTH + 0.02, self.SIMULATION_WIDTH, -self.SIMULATION_DEPTH)

        glEnd()

    def setup_points_rendering(self):
        """Setup GL state for smooth borderless points"""
        glEnable(GL_POINT_SMOOTH)  # Enable smooth points
        glEnable(GL_BLEND)         # Enable blending
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glDisable(GL_LIGHTING)     # Disable lighting
        glDisable(GL_DEPTH_TEST)   # Disable depth testing for particles
        
        # Set point parameters
        glPointSize(PARTICLE_SIZE)
        glHint(GL_POINT_SMOOTH_HINT, GL_NICEST)

        # Optional: Enable point sprites for better quality
        if GL_POINT_SPRITE in OpenGL.GL.__dict__:
            glEnable(GL_POINT_SPRITE)

    def draw_3d_dam(self):
        glBegin(GL_QUADS)
        glColor4f(1.0, 0.0, 1.0, 1)
        
        glVertex3f(-0.08, -self.SIMULATION_WIDTH, -self.SIMULATION_WIDTH)
        glVertex3f(0, -self.SIMULATION_WIDTH, -self.SIMULATION_WIDTH)
        glVertex3f(0, self.SIMULATION_WIDTH, -self.SIMULATION_WIDTH)
        glVertex3f(-0.08, self.SIMULATION_WIDTH, -self.SIMULATION_WIDTH)
        
        glVertex3f(-0.08, -self.SIMULATION_WIDTH, self.SIMULATION_WIDTH)
        glVertex3f(0, -self.SIMULATION_WIDTH, self.SIMULATION_WIDTH)
        glVertex3f(0, self.SIMULATION_WIDTH, self.SIMULATION_WIDTH)
        glVertex3f(-0.08, self.SIMULATION_WIDTH, self.SIMULATION_WIDTH)
        
        glColor4f(1.0, 0.5, 0.0, 0.7)
        glVertex3f(-0.08, 0, -self.SIMULATION_WIDTH)
        glVertex3f(-0.08, 0, self.SIMULATION_WIDTH)
        glVertex3f(-0.08, self.SIMULATION_WIDTH, self.SIMULATION_WIDTH)
        glVertex3f(-0.08, self.SIMULATION_WIDTH, -self.SIMULATION_WIDTH)
        glEnd()

    def color_velocity(self, particle):
        velocity = sqrt(particle.x_vel**2 + particle.y_vel**2)
        if velocity > 0.03:
            glColor4f(1.0, 0.0, 0.0, 0.8)
        elif velocity > 0.01:
            glColor4f(1.0, 1.0, 0.0, 0.8)
        else:
            glColor4f(0.0, 0.0, 1.0, 0.8)

    def color_press(self, particle):
        pressure = particle.press
        if pressure > 0.0002:
            glColor4f(.8, 0.0, 0.0, 0.8)
        elif pressure > 0.00005:
            glColor4f(.8, .8, 0.0, 0.8)
        else:
            glColor4f(0.0, 0.0, 0.8, 0.8)
    
    def draw_fps(self):
        fps = int(self.clock.get_fps())
        fps_text = self.font.render(f'FPS: {fps} | Frame Time: {(1000 / (fps + 1)):.4f} MEM: {process.memory_info().rss} Bytes', True, pygame.Color('white'))
        fps_data = pygame.image.tostring(fps_text, 'RGBA', True)
        glWindowPos2d(10, 10)  # Position at (10, 10)
        glDrawPixels(fps_text.get_width(), fps_text.get_height(), GL_RGBA, GL_UNSIGNED_BYTE, fps_data)

    def update(self):
        if self.frame == self.DAM_BREAK_FRAME:
            print("Breaking the dam")
            self.dam_built = False

        self.simulation_state = self.update_simulation(self.simulation_state, self.dam_built)
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()
        
        self.camera_x = self.camera_distance * np.sin(np.radians(self.rot_y))
        self.camera_z = self.camera_distance * np.cos(np.radians(self.rot_y))
        gluLookAt(
            self.camera_x, self.camera_distance * np.sin(np.radians(self.rot_x)), self.camera_z,
            0, 0, 0,
            0, 1, 0
        )

        if self.dam_built:
            self.draw_3d_dam()

        self.draw_3d_walls()

        glBegin(GL_POINTS)
        for particle in self.simulation_state:
            # change color based on the density of the particle
            self.color_func(particle)
            glVertex3f(
                particle.visual_x_pos,
                particle.visual_y_pos,
                particle.visual_z_pos
            )
        glEnd()

        self.frame += 1

        self.draw_fps()
        pygame.display.flip()
        self.clock.tick()

    def update_simulation(self, particles, dam):
        [p.update_state(dam) for p in particles]

        calculate_density(particles)

        [p.calculate_pressure() for p in particles]

        create_pressure(particles)
        calculate_viscosity(particles)
        
        return particles

    def start(self):
        while self.running:
            self.handle_input()
            self.update()
        pygame.quit()

if __name__ == "__main__":
    process = psutil.Process()
    animations = AnimatedScatter()
    animations.start()