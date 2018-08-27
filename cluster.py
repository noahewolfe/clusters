from amuse.units import nbody_system
from amuse.units import units
from amuse.ic.plummer import new_plummer_sphere
#from amuse.datamodel import Particles
from amuse.ic.salpeter import new_salpeter_mass_distribution

class Cluster:

	def __init__(self, num_stars, xoffset, yoffset, zoffset, repopulate_every_run=True):
		self.num_stars = num_stars
		self.xoffset = xoffset
		self.yoffset = yoffset
		self.zoffset = zoffset

		self.repopulate_every_run = repopulate_every_run
		if repopulate_every_run is False:
			# create nbody_converter thing?
 			convert_nbody = nbody_system.nbody_to_si(100.0 | units.MSun, 1 | units.parsec)
		
			# instantiate plummer sphere
			self.initial_plummer = new_plummer_sphere(self.num_stars, convert_nbody)
			self.initial_plummer.x += self.xoffset | units.lightyear
			self.initial_plummer.y += self.yoffset | units.lightyear
			self.initial_plummer.z += self.zoffset | units.lightyear

			# provide mass distribution for stars
			self.initial_plummer.mass = new_salpeter_mass_distribution(self.num_stars)

		self.plummer = None

	def populate(self):

		if self.repopulate_every_run:
			# create nbody_converter thing?
 			convert_nbody = nbody_system.nbody_to_si(100.0 | units.MSun, 1 | units.parsec)
		
			# instantiate plummer sphere
			self.plummer = new_plummer_sphere(self.num_stars, convert_nbody)
			self.plummer.x += self.xoffset | units.lightyear
			self.plummer.y += self.yoffset | units.lightyear
			self.plummer.z += self.zoffset | units.lightyear

			# provide mass distribution for stars
			self.plummer.mass = new_salpeter_mass_distribution(self.num_stars)
		else:
			self.plummer = self.initial_plummer

