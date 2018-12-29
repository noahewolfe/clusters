# clusters
## About
Hello! This project is currently under construction by Noah Wolfe. If you'd like to contact me, just email me at newolfe@ncsu.edu or noahwolfe1999@gmail.com.

The main goal of this project is to create a simulation using the [AMUSE framework](http://amusecode.org) in python, which simulates the collision of two open star clusters. This is to determine the types of complex structures which may emerge due to gravitational interactions on local interstellar scales. It is also to determine if the presence of stars, especially Blue Stragglers, greater than the 150 stellar mass limit (Crowther et al., 2010) could be due to collisions of stars during a cluster merger.

This repository was migrated from [this previous, somewhat messier, repository of mine](https://www.github.com/thezenth/cluster-collision)

## Running
This is currently a work in progress/mostly a reference for myself and my environment at the moment!

1. Write a new, unique simulation name in the "name" field in sim.yaml
2. Execute the following command from the "clusters" directory, assuming you have amuse installed in `~/amuse`: `mpiexec ~/amuse/amuse.sh main.py` (you need to execute it from the clusters directory, as it currently looks for sim.yaml in the working directory)

References

Crowther, Paul A., Olivier Schnurr, Raphael Hirschi, Norhasliza Yusof, Richard J. Parker, Simon P. Goodwin, and Hasan Abu Kassim. "The R136 Star Cluster Hosts Several Stars Whose Individual Masses Greatly Exceed the Accepted 150 M⊙ Stellar Mass Limit." Monthly Notices of the Royal Astronomical Society 408.2 (2010): 731-51. Web. 9 May 2016.
