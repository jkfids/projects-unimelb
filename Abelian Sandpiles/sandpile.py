import numpy as np
import math as m


class SandPile:
    """SandPile class
    """

    def __init__(self, width, height, threshold=4):
        """Initialize a sandpile with the specified width and height."""
        self.width = width
        self.height = height
        self.threshold = threshold
        self.grid = np.zeros((width, height), dtype=int)
        
        # Initializes a toppled grid which will be used to represent which
        # sites have been affected by an avalanche
        self.toppled_grid = np.zeros((width, height), dtype=int)
        
        # Initialize various properties such as histories for avalanche
        # properties, mass history and time
        self.time = 0;
        self.mass_history = [self.mass()]
        
        self.grain_loss = 0
        self.grain_loss_history = []
        self.grain_loss_total = 0
        
        self.topples_history = []
        self.area_toppled_history = []
        self.M_length_history = []
        self.E_length_history = []
        

    def drop_sand(self, n=1, site=None):
        """"Add `n` grains of sand to the grid.  Each grains of sand is added to
        a random site.
        
        This function also increments the time by 1 and update the internal
        `mass_history`.  Depending on how you want to code things, you may wish
        to also run the avalanche (alternatively, the avalanching might be
        executed elsewhere).

        Parameters
        ==========

        n: int

          The number of grains of sand of drop at this time step.  If left
          unspecified, defaults to 1.

        site:

          The site on which the grain(s) of sand should be dropped.  If `None`,
          a random site is used."""
        
        # Generates a random site when non is inputted
        if site==None:
            x = np.random.randint(low=0, high=self.width);
            y = np.random.randint(low=0, high=self.height);
            
        # Assumes that the site is given by 1x2 array where the elements are
        # the row and column number respectively
        else:
            x = site[0];
            y = site[1];
        
        # Increments the site by the number of sand grains inputted
        self.grid[x,y] = self.grid[x,y] + n;
        
        # Avalanche is run so that all unstable sites are stabilized
        self.avalanche();
        
        # Appends the length histories of the avalanche using our length functions
        # Also increments time and appends the mass history
        self.M_length_history.append(self.M_length_max(x,y,self.toppled_grid))
        self.E_length_history.append(self.E_length_max(x,y,self.toppled_grid))
        self.time = self.time + 1;
        self.mass_history.append(self.mass());

    def mass(self):
        """Return the mass of the grid using a double sum function"""
        return sum(sum(self.grid));
    
    def area_toppled(self):
        """Return the area toppled given the toppled_grid using a double sum"""
        return sum(sum(self.toppled_grid));
    
    def M_length(self, x1, y1, x2, y2):
        """Return the Manhatten distance between sites [x,y] and [i,j]"""
        return abs(x1-x2) + abs(y1-y2);
    
    def E_length(self, x1, y1, x2, y2):
        """Return the Euclidean distance between sites [x,y] and [i,j]"""
        return m.sqrt((abs(x1-x2)**2)+(abs(y1-y2))**2);
        
    def M_length_max(self, x, y, toppled_grid):
        """Return the M distance between the furthest affected site due to an 
        avalanche and the grain drop site that caused the avalanche"""
        
        length_max = 0;
        # Uses a double for loop to check the distance of each affected point
        # to the origin site and checks if it is the highest length
        for i in range(self.width):
            for j in range(self.height):
                if toppled_grid[i,j] == 1:
                    length = self.M_length(x,y,i,j);
                    if length_max < length:
                        length_max = length;
                        
        return length_max;
    
    def E_length_max(self, x, y, toppled_grid):
        """Return the E distance between the furthest affected site due to an 
        avalanche and the grain drop site that caused the avalanche"""
        
        length_max = 0;
        # Similiar to M_length_max()
        for i in range(self.width):
            for j in range(self.height):
                if toppled_grid[i,j] == 1:
                    length = self.E_length(x,y,i,j);
                    if length_max < length:
                        length_max = length;
                        
        return length_max;

    def topple(self, site):
        """Topple the specified site."""
        
        # Uses same site definition as in drop_sand()
        x = site[0];
        y = site[1];
        
        # Decreases the unstable site element by 4 and changes the respective
        # toppled_grid element to 1
        self.grid[x,y] = self.grid[x,y] - 4;
        self.toppled_grid[x,y] = 1;
        
        
        # Uses 4 if statements to increment each of the four adjacent sites
        # by 1 and adjusts the respective grid elemenent. However if the site
        # is not within the bounds of the grid then the grain loss count is
        # incremented by 1.
        if x < (self.width-1):
            self.grid[x+1,y] = self.grid[x+1,y] + 1;
            self.toppled_grid[x+1,y] = 1;
        else: self.grain_loss = self.grain_loss + 1;
        
        if x > 0:
            self.grid[x-1,y] = self.grid[x-1,y] + 1;
            self.toppled_grid[x-1,y] = 1;
        else: self.grain_loss = self.grain_loss + 1;
        
        if y < (self.height-1):
            self.grid[x,y+1] = self.grid[x,y+1] + 1;
            self.toppled_grid[x,y+1] = 1;
        else: self.grain_loss = self.grain_loss + 1;
        
        if y > 0:
            self.grid[x,y-1] = self.grid[x,y-1] + 1;
            self.toppled_grid[x,y-1] = 1;
        else: self.grain_loss = self.grain_loss + 1;
        

    def avalanche(self):
        """Run the avalanche causing all sites to topple and store the stats of
        the avalanche in the appropriate variables.
        """
        
        # Initializes avalanche properties that are recorded after the sandpile
        # is stabilized
        topples = 0;
        self.grain_loss=0;
        self.toppled_grid = np.zeros((self.width, self.height), dtype=int)
        
        # While there are any unstable sites, the function uses a double for
        # loop to scan through the entire grid and topple every unstable site 
        # Also increments the topples count each time a toppling occurs
        while np.any(self.grid >= 4):
            for i in range(self.width):
                for j in range(self.height):
                    if self.grid[i,j] >= 4:
                        self.topple([i,j]);
                        topples = topples + 1;
        
        # Appends various avalanche properties to the appropriate array histories
        self.topples_history.append(topples);
        self.grain_loss_history.append(self.grain_loss);
        self.grain_loss_total = self.grain_loss_total + self.grain_loss;
        self.area_toppled_history.append(self.area_toppled())
        
        