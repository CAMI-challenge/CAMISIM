#!/usr/env/python

with open('ecoli_aligned_length_ecdf','r') as length:
	aligned_bins = []
	for line in length:
		if not line.startswith('bin'):
			aligned_bins.append(line.split())

for length_bin in aligned_bins:
	length_bin[1] = float(length_bin[1])
	min_max = length_bin[0].split('-')
	length_bin[0] = (int(min_max[1]) - int(min_max[0]))/2 + int(min_max[0])

previous = 0.
mean = 0.
for length, running_sum in aligned_bins:
	mean += length * (running_sum - previous)
	previous = running_sum

print mean
