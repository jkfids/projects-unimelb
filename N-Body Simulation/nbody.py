# Template N-Body code by Rob Barone-Nugent
# Adapted from C code written by Chris Fluke and Matt Whiting
# December 2015: Third Year Astrophysics Labs

# The code provided here is a guide for those who are  
# unfamiliar with Python programming.  To learn more about 
# the different functions, statements, etc. introduced 
# here, you should have a look at some introductory Python
# documentation (e.g. https://wiki.python.org/moin/BeginnersGuide
# or ask your demonstrator....     

# Solving any problem with computers requries the
# following sequence of steps:

#     1) Determine exactly what it is you are
#            trying to solve.                     
#    2) Write the algorithm (eg. flow charts)  
#    3) Write the code.               
#    4) Run the code.                    
#    5) Test the code.               
#    6) Use the final code to get the answer.   
#                                                 
# It is best to get small sections of the code    
# working at a time, rather than trying to write  
# the entire program in one go -- this will make  
# finding errors a lot easier, as they are bound 
# to be there....                                

# This code provides the basic structure for an N-body program
# where the forces are calculated between each pair of particles 
# which provides accelerations and hence changes to the velocities
# and positions.  As written, the code will work for 2 bodies, and 
# in 2 dimensions (x, y) only.  You will need to make some changes
# to use the program for 3 bodies in 3 dimensions.                 

import numpy as np

# 'numpy' is a Python module, which provides a library of 
# numerical functions and structures for use in Python programs
# In order to     

m_unit = 2.e30        # Mass unit in kg
l_unit = 1.496e11    # Length unit in m
t_unit = 86400.        # Time unit in s , ie. second per day
G = 6.67e-11*m_unit*(t_unit**2)/(l_unit**3)        # The gravitational constant in rescaled units

# We define constants at the start of the program
# It is good programming practice to avoid using "magic numbers"  
# in your code, in case you decide later that you want to change 
# the value of one of a particular constant.    
# You WILL need to change the values of these constants

n_body = 4             # Number of bodies in the simulation
max_t = 50000            # Total simulation time
delta_t = 1         # Length of time step (in time units)
max_step = int(max_t/delta_t)    # Maximum number of time steps
show_step = 10    # Show every nth step

# We will use a Python structure called a 'class' to define the bodies
# The class will contain the mass, positions and velocities of every
# object
class body:
    def __init__(self, name, m, x, y, vx, vy, A):
        self.name = name    
        self.mass = m
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.A = A # Albedo
        
        self.orbits = None # Number of orbits around the Sun
        self.temp = None # Surface temperature of the planet
        self.phi = None # Current angle of orbit around the Sun
        self.phi_i = None # Angle of orbit before the Euler update
        self.phi_o = None # Original angle at time = 0
    # we will also add a function to the class that displays the 
    # body's parameters
    def show(self):
        print("Body: %s"%self.name)
        print("Mass = %6.2e mass units = %6.2e kg"%(self.mass, self.mass*m_unit))
        print("Position = (%6.2e,%6.2e) length units"%(self.x, self.y))
        print("Velocity = (%6.2e,%6.2e) velocity units"%(self.vx,self.vy))
        print("Kinetic Energy = %6.2e energy units"%(self.K_energy()))
        if self.temp != None:
            print("Surface Temperature =", int(self.temp), "K")
        if self.orbits != None:
            print("Number of Orbits:", int(self.orbits))
        print("")
    # Calculates the distance between self and another body
    def r(self, body):
        return np.sqrt((self.x - body.x)**2 + (self.y - body.y)**2)
    # Calculates the kinetic energy of the body
    def K_energy(self):
        return 0.5*self.mass*(self.vx**2 + self.vy**2)
    # Calculates the gravitational potential energy between self and another body
    def V_energy(self, body):
        return -G*self.mass*body.mass/self.r(body)
    # Calculates the surface temperature of a Planet
    def calc_temp(self, Sun):
        return 279*((1-self.A)**0.25)*(self.r(Sun)**-0.5)
        

# Now that we have defined the global variables and the body class
# we are now ready to run the simulation

# We will now define the initial conditions of the bodies and create
# them as an object (using the 'body' class) and add them to the
# 'bodies' array

# define the initial conditions of the Sun
name1 = 'Sun'
mass1 = 1
x1 = 50
y1 = 0.
vx1 = 0
vy1 = np.sqrt(G*50/x1)
A1 = None

# define the initial conditions of the Earth
name2 = 'Earth'
mass2 = 3.e-6
x2 = x1 + 1.
y2 = y1 + 0.
vx2 = vx1 + 0.
vy2 = vy1 + np.sqrt(G*mass1/(x2-x1))
A2 = 0.4

# define the initial conditions of Jupiter
name3 = 'Jupiter'
mass3 = 1/1047
x3 = x1 + 5.2
y3 = y1 + 0.
vx3 = vx1 + 0.
vy3 = vy1 + np.sqrt(G*mass1/(x3-x1))
A3 = 0.51

# define the initial conditions of the Moon
"""name4 = 'Moon'
mass4 = 0.0123 * mass2
x4 = x1 + x2 + 2.569e-3
y4 = y1 + y2 + 0.
vx4 = vx1 + vx2 + 0.
vy4 = vy1 + vy2 + np.sqrt(G*mass2/(x4-x1-x2))"""

# define the initial conditions of the black hole
name4 = 'Black Hole'
mass4 = 50
x4 = 0
y4 = 0
vx4 = 0
vy4 = 0
A4 = None

# Update the initial velocity of the Sun so the net linear momentum of the solar system = 0
vx1 -= ((vx2-vx1)*mass2 + (vx3-vx1)*mass3)/mass1
vy1 -= ((vy2-vy1)*mass2 + (vy3-vy1)*mass3)/mass1
# Does the same for the blackhole
vx4 -= vx1*mass1/mass4
vy4 -= vy1*mass1/mass4

# Define an array that will contain the body classes
# For now, the array is empty, but we will add bodies to it as needed
bodies = np.array([])

# Add the first body to the 'bodies' array
bodies = np.append(bodies, body(name1, mass1, x1, y1, vx1, vy1, A1))
bodies[0].temp = 5778 # Manually define the surface temperature of the Sun
# and the second
bodies = np.append(bodies, body(name2, mass2, x2, y2, vx2, vy2, A2))
# and the third and so on
bodies = np.append(bodies, body(name3, mass3, x3, y3, vx3, vy3, A3))
bodies = np.append(bodies, body(name4, mass4, x4, y4, vx4, vy4, A4))

# Returns the polar components of the vector between 2 Cartesian coordinates
def cart2pol(x1, y1, x2, y2):
    dx = x2 - x1
    dy = y2 - y1
    rho = np.sqrt(dx**2 + dy**2)
    phi = np.arctan2(dy, dx)
    return(rho, phi)

# Returns the total energies of the system
def system_energy(n_body=n_body, bodies=bodies):
    kinetic = 0.
    potential = 0.
    # Utilizes a for loop to sum up the kinetic energies and a double for loop
    # and an if statement to sum the potential energy between each unique pair
    for i in range(0, n_body):
        kinetic += bodies[i].K_energy()
        for j in range(0, n_body):
            if j > i:
                potential += bodies[i].V_energy(bodies[j])
    return (kinetic + potential, kinetic, potential)

# Initialize system/planetary properties and histories
temp_history = np.array([])
initial_energies = system_energy()
for i in range(1, n_body-1):
    bodies[i].orbits = 0.
    bodies[i].temp = bodies[i].calc_temp(bodies[0])
    bodies[i].phi = cart2pol(bodies[0].x, bodies[0].y, bodies[i].x, bodies[i].y)[1]
    bodies[i].phi_o = bodies[i].phi
    bodies[i].phi_i = bodies[i].phi
temp_history = np.append(temp_history, bodies[1].temp)

# The next part of this script is a loop -- it means that the 
# piece of indented code will be executed a number of times, 
# as the variable i takes on values from 1 to n_body. 
    
# Initialize and print system properties
print ("\033[4mInitial Properties\033[0m")
for i in range(0, n_body):
    bodies[i].show()    # Display name mass, position and velocity (as per the function 'show()' in the class 'body')
print("System:")
print("Total Kinetic = %6.3e energy units"%initial_energies[1])
print("Total Potential = %6.3e energy units"%initial_energies[2])
print("Total Energy = %6.3e energy units"%initial_energies[0])
print("")

# Open a file for each body that we will write the output to
# By using the open() function with the "w" (write) argument, a new
# file will be created and opened at the path of the first argument
f = open("outfile.csv", "w")
f.write(", ".join(["Xpos body%i, Ypos body%i"%(i+1, i+1) for i in range(n_body)])+", step, time,\n")

# We are now ready to run the n-body simulation
time = 0.
for step in range(0,max_step):
    # For each step, we want to record the acceleration/force between each pair of bodies
    # We create the two-dimensional array and populate it with zeros to begin
    # We want to calculate this array every step, so we create it inside the loop
    aFactor = np.zeros([n_body, n_body])

    # Work out all the separations and determine the acceleration factor
    for i in range(0,n_body):
        for j in range(0,n_body):
            # The "if" statements asks the computer to make a comparison
            # between the current values of i and j.  "!=" means not
            # equal to, as we do not want to calculate a self force.
            if i != j: 
                xsq = (bodies[i].x - bodies[j].x)**2.
                ysq = (bodies[i].y - bodies[j].y)**2.
                rsq = xsq + ysq
                factor = rsq**-1.5
                # update the acceleration factor array
                aFactor[i][j] = G*bodies[j].mass*factor
                aFactor[j][i] = G*bodies[i].mass*factor

    # Note that when we set up the class, we do not initialize the acceleration of
    # each body. So now we add a new feature to each of the bodies - the acceleration.
    for i in range(0,n_body):
           bodies[i].ax = 0.            # Set the accelerations to 0
           bodies[i].ay = 0.

    # And update the acceleration for each pair of bodies
    for i in range(0,n_body):
        for j in range(0,n_body):
            if i != j:
                # For each body, calculate the acceleration vector components
                bodies[i].ax -= aFactor[i][j] * (bodies[i].x - bodies[j].x)
                bodies[i].ay -= aFactor[i][j] * (bodies[i].y - bodies[j].y)

    # Save the initial phi values of the planets' orbit
    for i in range(0, n_body):
        if bodies[i].orbits != None:
            bodies[i].phi_i = cart2pol(bodies[0].x, bodies[0].y, bodies[i].x, bodies[i].y)[1]

    # We now update the position and velocity values for each body
    for i in range(0,n_body):
        # For each body, calculate the new velocity vector components
        bodies[i].vx += bodies[i].ax * delta_t;
        bodies[i].vy += bodies[i].ay * delta_t;
        # and the new position vector components
        bodies[i].x  += bodies[i].vx * delta_t;
        bodies[i].y  += bodies[i].vy * delta_t;
    
    # Update the phi component of the planets' orbital path then checks whether
    # the body has completed an orbit
    # See labnotes for more details
    for i in range(0,n_body):
        if bodies[i].orbits != None:
            bodies[i].phi = cart2pol(bodies[0].x, bodies[0].y, bodies[i].x, bodies[i].y)[1]
            if bodies[i].phi_o < bodies[i].phi and bodies[i].phi_o > bodies[i].phi_i:
                bodies[i].orbits += 1
    
    # Update the planetary temperatures and append it to the temperature history
    for i in range(0,n_body):
        if bodies[i].A != None:
            bodies[i].temp = bodies[i].calc_temp(bodies[0])
    temp_history = np.append(temp_history, bodies[1].temp)
            
    # Update the system energies
    energies = system_energy()
    
    # We don't want to write out data every step, just every show_step
    # times -- EXCEL can't handle more than about 4000 values.         
    if (step%show_step)==0:
        for i in range(0,n_body):
            ## Write out the timestep, position and velocity data for each particle
            f.write("%0.10f,%0.10f,"%(bodies[i].x, bodies[i].y))
        f.write("%5d,%6.4f\n"%(step, time))

    # update the current time
    time+=delta_t

# Once the loop is finished, close the file so we can read it.
f.close()

# plot
# I also added legends and markers for the final positions of the bodies.
import matplotlib.pyplot as plt
plt.figure(figsize=(6,6), dpi=70)
plt.axis('equal') # Makes the scale of the two axis equal, so circles appear as circles
f = open("outfile.csv"); k = f.readlines(); f.close()
x_body1 = [float(i.split(',')[0]) for i in k[1:]]
y_body1 = [float(i.split(',')[1]) for i in k[1:]]
x_body2 = [float(i.split(',')[2]) for i in k[1:]]
y_body2 = [float(i.split(',')[3]) for i in k[1:]]
x_body3 = [float(i.split(',')[4]) for i in k[1:]]
y_body3 = [float(i.split(',')[5]) for i in k[1:]]
x_body4 = [float(i.split(',')[6]) for i in k[1:]]
y_body4 = [float(i.split(',')[7]) for i in k[1:]]
plt.plot(x_body1, y_body1, c='y')
plt.plot((x_body1[len(x_body1)-1]),(y_body1[len(y_body1)-1]), marker='o', c='y', label='Sun')
plt.plot(x_body2, y_body2)
plt.plot((x_body2[len(x_body2)-1]),(y_body2[len(y_body2)-1]), marker='o', c='#1f77b4', label='Earth')
plt.plot(x_body3, y_body3)
plt.plot((x_body3[len(x_body3)-1]),(y_body3[len(y_body2)-3]), marker='o', c='#ff7f0e', label='Jupiter')
plt.plot(x_body4, y_body4, c='black')
plt.plot((x_body4[len(x_body4)-1]),(y_body4[len(y_body4)-1]), marker='o', c='black', label='Black Hole')
plt.title('Black Hole vs Solar System - 50000 days (Δt = 1)')
plt.legend()
plt.show()

# Plot the surface temperature of Earth over the simulation
t_history = delta_t*np.array([i for i in range(len(temp_history))])
plt.plot(t_history, temp_history)
plt.title("Surface Temperature of Earth over Time")
plt.xlabel("Time (day)")
plt.ylabel("Temperature (K)")

# Prints properties of the system after simulation
print ("\033[4mFinal Properties\033[0m")
for i in range(0, n_body):
    bodies[i].show()    # Display name mass, position and velocity (as per the function 'show()' in the class 'body')
print("System:")
print("Total Kinetic = %6.3e energy units"%energies[1])
print("Total Potential = %6.3e energy units"%energies[2])
print("Total Energy = %6.3e energy units"%energies[0])