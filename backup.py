import os

outfiles = [i for i in os.listdir('.') if i.endswith('.csv')]
ids = []
real_ids = list(range(539))

for outfile in outfiles:
	job_id = (outfile.split('_')[1]).split('.')[0]
	myjob_id = int(job_id)
	ids.append(myjob_id)

for jobid in real_ids:
	if jobid not in ids:
		print(jobid)



