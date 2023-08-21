import os
import re

d = open('jobs.txt','r')
t = d.read()
jobs = [int(i) for i in t.split('\n')]
d.close()


outfiles = [i for i in os.listdir('.') if i.endswith('.out')]

for outfile in outfiles:
	d = open(outfile,'r')
	t = d.read()
	d.close()
	success = re.findall('(?<=FINISHED BATCH )[0-9]+',t)

	if len(success) == 1:
		batch_id = int(success[0])
		jobs.remove(batch_id)


d = open('jobs.txt','w')
d.write('\n'.join([str(i) for i in jobs]))
d.close()
