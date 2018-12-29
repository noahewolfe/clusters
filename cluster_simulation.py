import os
import math

import h5py
from amuse.units import nbody_system
from amuse.units import units
from amuse.datamodel import Particles
from amuse.community.hermite0.interface import Hermite
from amuse.community.sse.interface import SSE

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
import yaml

OUTPUT_PATH = "/home/noah/clusters_output/"
if os.path.exists(OUTPUT_PATH) is False:
	os.mkdir(OUTPUT_PATH)

"""def get_color_based_on_stellar_type(istp):
    color = get_distinct(5)
    if istp.value_in(units.stellar_type) < 2:
        c = color[0]
    elif istp.value_in(units.stellar_type) < 6:
        c = color[3]
    else:
        c = color[4]
return c"""

def plot_data(plot_name, plot_file_name, particles, zoom=0.0, bound=20.0, mass_u=units.MSun, dist_u=units.lightyear):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d') #make it a 3d plot

	xdata = []
	ydata = []
	zdata = []
	massdata = []

	for p in particles:
		xdata.append(p.position.x.value_in(dist_u))
		ydata.append(p.position.y.value_in(dist_u))
		zdata.append(p.position.z.value_in(dist_u))
		massdata.append(p.mass.value_in(mass_u))

	#print(xdata)
	#print(ydata)
	#print(zdata)

	ax.scatter(xdata, ydata, zdata, s=massdata)

	ax.set_xlabel(dist_u)
	ax.set_ylabel(dist_u)
	ax.set_zlabel(dist_u)

	ax.set_xlim([-bound, bound]) #set this bound as graph limit
	ax.set_ylim([-bound, bound])
	ax.set_zlim([-bound, bound])

	fig.suptitle(plot_name)

	plt.savefig(plot_file_name)
	#plt.show()
	plt.close()


def write_csv(file_name, particles, time, mass_u=units.kg, dist_u=units.lightyear):
	f = open(file_name, "w")
	f.write("TIME=" + str(time) + "\n")
	f.write("id,mass,x,y,z\n")

	for i in range(0, len(particles)):
		p = particles[i]

		ident = str(i)
		mass = str(p.mass.value_in(mass_u))
		xpos = str(p.position.x.value_in(dist_u))
		ypos = str(p.position.y.value_in(dist_u))
		zpos = str(p.position.z.value_in(dist_u))

		f.write(
			ident + "," +
			mass  + "," +
			xpos  + "," +
			ypos  + "," +
			zpos  + "\n"
		)

	f.close()

#def write_collision_csv(file_name, )


class ClusterSimulation:
	def __init__(self, name, timestep, runtime, num_runs, clusters):
		
		print("Initializing cluster simulation %s" % (name) )
		
		self.name = name
		self.num_runs = num_runs
		self.timestep = timestep
		self.runtime = runtime
		self.clusters = clusters

		max_iterations = math.ceil(timestep / runtime) # ceil-ed max number of iterations, to determine 0-padding in hdf5 group names
		self.iteration_name_length = len( str(max_iterations) )

		# initialize file system stuff
		# main folder for all runs
		print("		Creating main folder for run output...")
		self.folder = OUTPUT_PATH + self.name
		os.mkdir(self.folder)
		print("		Saving configuration file...")
		clusters_saved = [] # will contain objects with only the cluster info that we want... not backgronud amuse stuff
		for c in self.clusters:
			clusters_saved.append({
				"num_stars"            : c.num_stars,
				"xoffset"              : c.xoffset,
				"yoffset"              : c.yoffset,
				"zoffset"              : c.zoffset,
				"repopulate_every_run" : c.repopulate_every_run
			})

		### Save configuration ###
		data = dict(
			name = self.name,
			timestep = self.timestep,
			runtime = self.runtime,
			num_runs = self.num_runs,
			clusters = clusters_saved
		)
		with open(self.folder + "/config.yaml", "w") as outfile:
			yaml.dump(data, outfile, default_flow_style=False)
		print("		Done!")

		### Setup HDF5 File Object ###
		self.hdf5_file = h5py.File(self.folder + "/" + self.name + ".hdf5", "w")

		print("Done initializing simulation!")

	def write_to_hdf5_file(self, run_group, iteration, time, particles, mass_u=units.kg, dist_u=units.lightyear):
		### Create iteration group ###
		## Create group name 
		iter_group_name = str(iteration)
		add_zeroes = self.iteration_name_length - len(iter_group_name)
		for z in range(add_zeroes):
			iter_group_name = "0" + iter_group_name

		## Create group
		iter_group = run_group.create_group(iter_group_name)
		iter_group.attrs.create("time", time)

		key_list = []

		x_list = []
		y_list = []
		z_list = []

		#vx_list = []
		#vy_list = []
		#vz_list = []

		mass_list = []
		#radius_list = []

		### Add datasets for each particle
		for p in particles:
			key_list.append(p.key)
			x_list.append(p.position.x.value_in(dist_u))
			y_list.append(p.position.y.value_in(dist_u))
			z_list.append(p.position.z.value_in(dist_u))
			#vx_list.append(p.vx)
			#vy_list.append(p.vy)
			#vz_list.append(p.vz)
			mass_list.append(p.mass.value_in(mass_u))
			#radius_list.append(p.radius)	

		keys = np.array(key_list, dtype="uint64")

		x = np.array(x_list, dtype="float")
		y = np.array(y_list, dtype="float")
		z = np.array(z_list, dtype="float")

		#vx = np.array(vx_list)
		#vy = np.array(vy_list)
		#vz = np.array(vz_list)

		masses = np.array(mass_list)
		#radii = np.array(radius_list)

		all_particles = np.vstack( (keys, x, y, z, masses) )

		iter_group.create_dataset("particles", data=all_particles)

	def run(self, run_num):
		print("Setting up run %d..." % (run_num))

		# Create group for run
		run_group = self.hdf5_file.create_group(str(run_num))

		# create nbody converter thing?
	 	convert_nbody = nbody_system.nbody_to_si(100.0 | units.MSun, 1 | units.parsec)

		# initialize particle datamodel
		stars = Particles()
		
		print("		Initializing and populating clusters...")
		# populate/initialize clusters, add to particle model
		for c in self.clusters:
			c.populate()
			stars.add_particles(c.plummer)
		print("		Done!")

		# initialize codes
		# initialize hermite
		print("		Initializing hermite...")
		hermite_code = Hermite(convert_nbody)
		hermite_code.particles.add_particles(stars)
		detect_coll = hermite_code.stopping_conditions.collision_detection
		detect_coll.enable()
		print("		Done!")

		print("		Initializing SSE...")
		sse_code = SSE()
		sse_code.particles.add_particles(stars)

		print("Done!\n")

		# actually run the simulation!
		print("===== Run #%d =====" % (run_num))

		t = 0 # time
		i = 1 # iteration num
		c = 1 # collision num
		while(t < self.runtime):
			print("		Time (Myr): %d" % (t))

			print("			Simulating...")

			# evolve model
			hermite_code.evolve_model(t | units.Myr)
			sse_code.evolve_model(t | units.Myr)

			hermite_code.particles.copy_values_of_attribute_to("position", sse_code.particles)
			sse_code.particles.synchronize_to(hermite_code.particles)

			if detect_coll.is_set():
				print("Detected a collision!!")
				print(detect_coll.particles(0))
				coll_f = open(sub_folder + "/collision-" + str(c) + "_time" + str(t) + ".txt")
				coll_f.write( detect_coll.particles(0) )
				coll_f.close()
				c += 1
				# somehow put a log
				# also... eventually, "pause" the rest of the simulation and simulate the collision somehow

			#sse_code.particles.copy_values_of_attribute_to("mass", hermite_code.particles)
			#sse_code.particles.copy_values_of_attribute_to("stellar_type", hermite_code.particles)
			#sse_code.particles.copy_values_of_attribute_to("age", hermite_code.particles)
			#sse_code.particles.copy_values_of_attribute_to("luminosity", hermite_code.particles)
			#sse_code.particles.copy_values_of_attribute_to("temperature", hermite_code.particles)
			#sse_code.particles.copy_values_of_attribute_to("radius", hermite_code.particles)

			hermite_code.particles.mass = sse_code.particles.mass

			print("			Done.")

			# output to csv
			print("			Outputting to HDF5...")
			#file_name = sub_folder + "/data-" + str(i) + "_time-" + str(t) + ".txt"
			self.write_to_hdf5_file(run_group, i, t, stars)
			print("			Done.")

			#print("			Creating plot...")
			#plot_name = "System at " + str(t) + " Myr"
			#plot_path = sub_folder + "/plot-" + str(i) + "_time-" + str(t) + ".png"
			#plot_data(plot_name, plot_path, hermite_code.particles)
			#print("			Done.")

			t += self.timestep
			i += 1

		print("Run complete.\n")


	def do_simulation(self):
		for r in range(1, self.num_runs + 1):
			self.run(r)
		self.hdf5_file.close()