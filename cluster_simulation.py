import os

from amuse.units import nbody_system
from amuse.units import units
from amuse.datamodel import Particles
from amuse.community.hermite0.interface import Hermite
from amuse.community.sse.interface import SSE

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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


class ClusterSimulation:
	def __init__(self, name, timestep, runtime, num_runs, clusters):
		print("Initializing cluster simulation %s" % (name) )
		self.name = name
		self.num_runs = num_runs
		self.timestep = timestep
		self.runtime = runtime
		self.clusters = clusters

		# initialize file system stuff
		# main folder for all runs
		print("		Creating main folder for run output...")
		self.folder = OUTPUT_PATH + self.name
		os.mkdir(self.folder)
		print("		Done!")

		print("Done initializing simulation!")

	def run(self, run_num):
		print("Setting up run %d..." % (run_num))

		# create sub-folder for the run
		print("		Creating folder for results...")
		sub_folder = self.folder + "/run_" + str(run_num)
		os.mkdir(sub_folder)
		print("		Done!")

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

		t = 0
		while(t < self.runtime):
			print("		Time (Myr): %d" % (t))

			# evolve model
			hermite_code.evolve_model(t | units.Myr)
			sse_code.evolve_model(t | units.Myr)

			hermite_code.particles.copy_values_of_attribute_to("position", sse_code.particles)
			sse_code.particles.synchronize_to(hermite_code.particles)

			if detect_coll.is_set():
				print("Detected a collision!!")
				print(detect_coll.particles(0))
				# somehow put a log
				# also... eventually, "pause" teh rest of the simulation and simulate the collision somehow

			#sse_code.particles.copy_values_of_attribute_to("mass", hermite_code.particles)
			#sse_code.particles.copy_values_of_attribute_to("stellar_type", hermite_code.particles)
			#sse_code.particles.copy_values_of_attribute_to("age", hermite_code.particles)
			#sse_code.particles.copy_values_of_attribute_to("luminosity", hermite_code.particles)
			#sse_code.particles.copy_values_of_attribute_to("temperature", hermite_code.particles)
			#sse_code.particles.copy_values_of_attribute_to("radius", hermite_code.particles)


			hermite_code.particles.mass = sse_code.particles.mass

			# output to csv
			print("			Outputting to CSV...")
			file_name = sub_folder + "/data_time-" + str(t) + ".txt"
			write_csv(file_name, stars, t)
			print("			Done.")

			print("			Creating plot...")
			plot_name = "System at " + str(t) + " Myr"
			plot_path = sub_folder + "/plot_time-" + str(t) + ".png"
			plot_data(plot_name, plot_path, hermite_code.particles)
			print("			Done.")

			t += self.timestep

		print("Run complete.\n")


	def do_simulation(self):
		for r in range(1, self.num_runs + 1):
			self.run(r)