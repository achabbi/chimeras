import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--batches', type=int, required=True, help='Number of batches')
args = parser.parse_args()

file = open('jobs.txt', 'w')

for idx in range(args.batches - 1):
	file.write(str(idx))
	file.write('\n')

file.write(str(args.batches - 1))

file.close()

