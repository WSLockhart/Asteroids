'''
Current Problems / Goals:

* what are the global variables?

- restructure using after()
- use move() instead of erasing the whole canvas
- objects still getting stuck behind the walls?
- enable elastic scattering
- test grad(U) method of getting the acceleration
- inputs for gaussians
- trace option
- autozoom?
- units
- canvas outline (aesthetic)
- ability to save settings :o
--------------------------------------------------------------------------'''

import numpy as np      #math package
import tkinter as tk    #graphics package

#----------------------------------------------------------------------
# a body is characterized by its mass, radius, position vector, and velocity vector
class body:
    def __init__(self, m, r, x, v):
        self.mass = float(m)
        self.radius = float(r)
        self.position = np.array(x)
        self.velocity = np.array(v)
        self.momentum = m*self.velocity
        # keep track of my acceleration so I can return it at the end of each increment of the sim
        # see simulation.run() for details on this subtle but important issue
        self.acceleration = np.zeros(self.velocity.size)

    def set_position(self, vector):
        self.position = np.array(vector)

    def set_velocity(self, vector):
        self.velocity = np.array(vector)

    def set_momentum(self, vector):
        self.momentum = np.array(vector)

    def set_acceleration(self, vector):
        self.acceleration = np.array(vector)

    def edit_coord(self, index, value):
        self.position[index] = value


'''------------------------------------------------------------------------'''

# # a window from which you can launch different types of mainframes!
# class masterframe:
#     def __init__(self, master):
#         frame = tk.Frame(master, bg='gray70')
#         frame.pack(fill=tk.BOTH, expand=1)
#         btn_frame = tk.Frame(frame, bd=5, relief=tk.RAISED)
#         btn_frame.pack(padx=10,pady=10)
#         tk.Button(btn_frame, text="classical simulation", command=self.launch_classical).pack(side=tk.LEFT)
#         tk.Button(btn_frame, text="quantum simulation", state='disabled',command=self.launch_quantum).pack(side=tk.LEFT)
#
#
#     def launch_classical(event=None):
#         new_root = tk.Tk()
#         new_root.title('CLASSICAL')
#         new_mf = mainframe(new_root)
#         new_root.mainloop()
#
#     def launch_quantum(event=None):
#         new_root = tk.Tk()
#         new_root.title('QUANTUM')
#         new_mf = mainframe(new_root)
#         new_root.mainloop()


# ---------------------------------------------------------------------- #
#                               GRAPHICS
# ---------------------------------------------------------------------- #
class mainframe:

    def __init__(self, master):
        base = tk.Frame(master, bg='gray70')
        base.pack(fill=tk.BOTH, expand=1)

        canvas_frame = tk.Frame(base, bg='gray70')
        canvas_frame.pack(side=tk.LEFT)
        self.scale_factor = 300   # determines the zoom
        self.canvas = tk.Canvas(canvas_frame, width=600, height=600, bd=5, relief=tk.RAISED)
        self.canvas.pack(padx=10, pady=15)
        btn_frame = tk.Frame(canvas_frame, bd=5, relief=tk.RAISED)
        btn_frame.pack(pady=10)
        tk.Button(btn_frame, text="load simulation", font="-weight bold", command=self.launch_new_sim).pack(side=tk.LEFT)
        tk.Button(btn_frame, text="run", command=self.run).pack(side=tk.LEFT)
        tk.Button(btn_frame, text="stop", command=self.stop).pack(side=tk.LEFT)
        tk.Button(btn_frame, text="zoom out", command=self.zoom_out).pack(side=tk.LEFT)
        tk.Button(btn_frame, text="zoom in", command=self.zoom_in).pack(side=tk.LEFT)
        tk.Button(btn_frame, text="exit", command=self.quit).pack(side=tk.LEFT)

        input_frame = tk.Frame(base, bd=5, relief=tk.RAISED)
        input_frame.pack(side=tk.LEFT, padx=10, pady=15)
        global D
        D = tk.IntVar()
        D.set('2')
        entry_frame_d = tk.Frame(input_frame)
        entry_frame_d.pack(pady=4)
        tk.Entry(entry_frame_d, width=6, textvariable=D, state='disabled').pack(side=tk.RIGHT)
        tk.Label(entry_frame_d, text='dimensions:').pack(side=tk.LEFT)
        global n
        n = tk.IntVar()
        n.set('20')
        entry_frame_n = tk.Frame(input_frame)
        entry_frame_n.pack(pady=4)
        tk.Entry(entry_frame_n, width=6, textvariable=n).pack(side=tk.RIGHT)
        tk.Label(entry_frame_n, text='bodies:').pack(side=tk.LEFT)
        # global G
        # G = tk.StringVar()
        global G_entry
        entry_frame_g = tk.Frame(input_frame)
        entry_frame_g.pack(pady=4)
        G_entry = tk.Entry(entry_frame_g, width=6)
        G_entry.insert(0,'0.001')
        G_entry.pack(side=tk.RIGHT)
        tk.Label(entry_frame_g, text='coupling:').pack(side=tk.LEFT)
        global sticky
        sticky = tk.IntVar()
        tk.Checkbutton(input_frame, text=' elastic', variable=sticky).pack(pady=4)
        global walled
        walled = tk.IntVar()
        walled.set('0')
        tk.Checkbutton(input_frame, text=' walls', variable=walled).pack(pady=4)
    # ----


    def launch_new_sim(event=None):
        print ('new simulation.')
        global current_sim  #needs to be global so that run() and stop() can access it

        # collect input from the input boxes:
        try: nb = n.get()
        except ValueError: print('# bodies must be a positive integer')
        try: g = float(G_entry.get())
        except ValueError: print('interaction strength must be a valid number')
        if sticky.get() == 1: elast = True
        else: elast = False
        if walled.get() == 1: walls = True
        else: walls = False

        current_sim = simulation(2,g,nb,elast,walls) # dim, grav, #bodies, elastic, walls
    # -----


    def run(event=None):
        current_sim.run()

    def stop(event=None):
        current_sim.stop()

    # Define a function to close the window.
    def quit(event=None):
        root.destroy()

    def zoom_out(self):
        self.scale_factor /= 2
        current_sim.draw_asteroids(current_sim.asteroids)


    def zoom_in(self):
        self.scale_factor *= 2
        current_sim.draw_asteroids(current_sim.asteroids)


    def draw_circle(self,x_coord,y_coord,rad, color):
        # convert from internal coords (-1 to 1) to canvas coords (0 to 600)
        x = 300 + self.scale_factor*x_coord
        y = 300 + self.scale_factor*y_coord
        r = rad*self.scale_factor
        self.canvas.create_oval(x-r,y-r,x+r,y+r, fill=color)

    def clear_canvas(self):
        self.canvas.delete('all')


'''------------------------------------------------------------------------'''


# simulation of N bodies in a D-dimensional box of sidelength 2
class simulation():

    def __init__(self, dimension, gravity_strength, num_bodies, elastic, walls):
        self.load_globals(dimension, gravity_strength, num_bodies, elastic, walls)
        #create randomly placed asteroids, mass range 10 <= m <= 100
        self.asteroids = self.load_asteroids(num_bodies, 10.0 , 100.0)
        self.draw_asteroids(self.asteroids)
        global sim_running
        sim_running = False


    # ----------------------------------------------------------------------
    def run(self):
        global sim_running
        sim_running = True
        dt = float(0.005)   #NOTE: dt determines the speed & accuracy of simulation
        print ('running sim...')
        while sim_running:
            for body in self.asteroids:
                # total_momentum = np.zeros(D)
                #1) Find the acceleration: a = F/m = -GMgrad(U). #NOTE: may want to use grad(U) instead
                acceleration = self.get_force(body) / body.mass
                #NOTE!!! CANNOT update the positions/velocities immediately or things will get off.
                # Object #1 will have a 'head start' that grows over time until when two bodies merge
                # they have a wildly incorrect momentum! Ah!
                body.set_acceleration(acceleration)

                #2) If you hit a wall, reverse direction
                if walls:
                    for i in range(body.position.size):
                        if abs(body.position[i]) > 1:
                            body.velocity[i] *= -1
                            body.edit_coord(i,round(body.position[i]))  #this is to avoid things getting stuck behind the walls

                #3) Check for collisions
                for otherbody in self.asteroids:
                    if np.linalg.norm(body.position - otherbody.position) < (body.radius + otherbody.radius)*0.8 and otherbody != body:
                        # if elastic:
                        #     # Momentum is conserved: m1v1 + m2v2 = MV
                        #     # Energy is conserved: (1/2)m1v1^2 + (1/2)m2v2^2 = (1/2)MV^2
                        # else:
                        newbody = self.merge_bodies(body,otherbody)
                        self.asteroids.remove(body)
                        self.asteroids.remove(otherbody)
                        self.asteroids.append(newbody)
                        break  #important! took forever to debug without this. Otherwise it thinks there's another merger
                        # between the new body and one of the old ones, even though the old one is gone!
            #| end for loop

            # Now update all the positions & velocities
            for body in self.asteroids:
                body.set_velocity( body.velocity + (body.acceleration*dt) )
                body.set_position( body.position + (body.velocity*dt) )
                body.set_momentum(body.mass*body.velocity)

            self.draw_asteroids(self.asteroids) #update the canvas
            root.update()  #refresh the tkinter graphics env

            # for ast in self.asteroids: total_momentum += ast.momentum
        #| end while loop


    def stop(self):
        global sim_running
        sim_running = False
    # ----------------------------------------------------------------------


    def load_globals(self, dimensions, gravity_strength, num_bodies, stic, bounded):
        #number of spatial dimensions
        global D
        D = int(dimensions)
        # strength of gravity
        global G
        G = float(gravity_strength)
        # N-body problem
        global N
        N = int(num_bodies)
        global elastic
        elastic = stic
        global walls
        walls = bounded
        #NOTE: asteroids is a variable of simulation right now, not a global
        # should all of these be self.variable instead?


    def load_asteroids(self, N, min_mass, max_mass):
        asteroid_list = []
        for i in range(N):
            new_asteroid = self.generate_random_body(min_mass,max_mass)
            # new_asteroid.set_velocity(np.zeros(D))
            asteroid_list.append(new_asteroid)
        return asteroid_list


    def draw_asteroids(self,bodies):
        mf.clear_canvas()
        for body in bodies:
            # small masses are darker shades
            colorcode = 'gray' + str(min(int(body.mass/2),99))
            mf.draw_circle(body.position[0],body.position[1],body.radius, colorcode)


    def generate_random_body(self,min_mass, max_mass):
        mass = np.random.uniform(min_mass,max_mass)
        radius = (mass**(1/D))/300  #volume proportional to mass. 600 is the mf's initial scale_factor
        position_vector = self.generate_gaussian_vector(D,0,0.4)  #NOTE: position distribution determined here
        velocity_vector = self.generate_random_vector(D,-1,1)
        return body(mass,radius,position_vector,velocity_vector)


    # return a new body. sum the mass, sum the velocities, average the positions
    def merge_bodies(self,body1,body2):
        new_mass = body1.mass + body2.mass
        new_radius = new_mass**(1/D)/300  #volume proportional to mass
        #NOTE: Another annoying bug! Don't just average positions or you'll introuce error
        # Must use the Center of Mass: COM = (m1r1 + m2r2)/(m1+m2)
        new_position = (body1.mass*body1.position + body2.mass*body2.position) / new_mass
        # Conservation of Momentum: m1v1 + m2v2 = (m1+m2)V, so V = (m1v1 + m2v2)/(m1+m2)
        new_velocity = (body1.mass*body1.velocity + body2.mass*body2.velocity)/(new_mass)

        return body(new_mass,new_radius,new_position,new_velocity)


    def generate_gaussian_vector(self,n,mean,sd):
        v = np.empty(n)
        for i in range(n):
            v[i] = np.random.normal(mean,sd)
        return v


    def generate_random_vector(self,n,_min,_max):
        v = np.empty(n)
        for i in range(n):
            v[i] = np.random.uniform(_min,_max)
        return v


    def get_distance(self,r1,r2):
        return np.linalg.norm(r1 - r2)


    def get_force(self,me):
        F = np.zeros(D)
        for otherbody in self.asteroids:
            if otherbody != me:
                # this r is a vector
                r = otherbody.position - me.position    #vector points from me to him
                # F = GMmr / |r|^3
                F += (G * me.mass * otherbody.mass) * r / (np.linalg.norm(r))**3
        return F


    # potential = sum[Gm/r]
    def get_potential(self,coordinates):
        U = float(0)
        for body in self.asteroids:
            # this r is a scalar
            r = self.get_distance(coordinates,body.position)
            if r > 0:  #avoids inifinte potential of self
                U += ((-1) * G * body.mass) / r
        return U


    # the field is the gradient of the potential: (d/dx,d/dy,d/dz)
    def gradient(self,position_vector):
        coordinates = np.array(position_vector)
        # to approximate the gradient, go out a little ways in each direction and take the slope, dU/dx
        grad = np.empty(D)
        dx = 0.01   #NOTE: contributes to the accuracy of the simulation. Right now it's set to 1/100 the size of the box
        new_coordinates = coordinates
        for i in range(2):
            new_coordinates[i] += dx
            grad[i] = (self.get_potential(coordinates) - self.get_potential(new_coordinates)) / dx
            new_coordinates[i] -= dx
        return grad


# -------------------  End of class simulation  ------------------- #


'''------------------------------------------------------------------------'''


# def main():
root = tk.Tk()
root.title('PHYSICS SIMULATOR')
mf = mainframe(root)
root.mainloop()


# #----------------------------------------------------------------------
# if __name__ == "__main__":
#     main()
