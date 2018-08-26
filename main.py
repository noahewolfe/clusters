from cluster import Cluster
from cluster_simulation import ClusterSimulation

import yaml

with open("./sim.yaml") as stream:
	try:
		config = yaml.load(stream)
		
		cluster_objs = []
		for c_desc in config["clusters"]:
			cluster_objs.append(
				Cluster( 
					c_desc["num_stars"],
					c_desc["xoffset"],
					c_desc["yoffset"],
					c_desc["zoffset"] 
				)
			)

		sim = ClusterSimulation(
			config["name"],
			config["timestep"],
			config["runtime"],
			config["num_runs"], 
			cluster_objs
		)

		sim.do_simulation()

		#print(config)
	except yaml.YAMLError as exc:
		print(exc)

"""
name = raw_input("Name: ")

num_runs = int( input("# of runs: ") )
while num_runs <= 0:
	print("# of runs must be greater than 0!")
	num_runs = int( raw_input("# of runs: ") )

timestep = float( raw_input("Timestep (Myr.): ") )
total_runtime = float( raw_input("Total Runtime (Myr.): ") )

num_clusters = int( raw_input("# of clusters: ") )
while num_clusters <= 0:
	print("# of clusters must be greater than 0!")
	num_clusters = int( raw_input("# of clusters: ") )

clusters = []

for c in range(1, num_clusters + 1):
	print("Parameters for cluster %d" % (c))

	num_stars = int( raw_input("# of stars: ") )
	xoff = int( raw_input("x offset (ly.): ") )
	yoff = int( raw_input("y offset (ly.): ") )
	zoff = int( raw_input("z offset (ly.): ") )

	clusters.append( Cluster(num_stars, xoff, yoff, zoff) )

sim = ClusterSimulation(name, timestep, total_runtime, num_runs, clusters)
sim.do_simulation()
"""