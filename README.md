# chimeras

Scripts for generating AAV protein chimeras using the SCHEMA recombination algorithm.
* `create_crossover_points.py` makes the crossover points where the sequences are swapped
* `chimeras.py` recombines the sequences and calculates the SCHEMA disruptions and hamming distances
* `backup.py`, `cleanup.py`, and `jobs.py` keep track of jobs on the cluster
