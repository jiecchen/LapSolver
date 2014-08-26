#!/usr/bin/python
import random
import sys

size = int(sys.argv[1])

nums = range(size)
random.shuffle(nums)
iNums = range(size)

for i in xrange(size):
    iNums[nums[i]] = i

with open('test.in', 'w') as f:
    f.write('{}\n'.format(size))
    for i in xrange(size):
        f.write('{} {}\n'.format(iNums[i]+1, i+1))

with open('test.out', 'w') as f:
    for i in nums:
        f.write('{}\n'.format(i+1))

