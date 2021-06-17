from os import listdir
from os.path import isfile, join, isdir
from subprocess import check_output

from functools import partial
from multiprocessing.dummy import Pool
from subprocess import call

import logging

Gs = ['g150.4.col', 'g150.5.col', 'g170.3.col', 'g200.2.col',
      'g200.3.col', 'g300.2.col', 'g350.2.col', 'g400.1.col']

# Gx TEST
DIRTEST = './data/sparse/'
COMMANDS = []

#for T in [11, 12, 13, 14, 5]:
#	for dt in Gs:
#    		TEST = join(DIRTEST, dt)
#
#    		LOG = join('./tmp/{}.test.{}.log'.format(dt, T))
#    		CMD = './bin/mssRank ' + TEST + ' {} > '.format(T) + LOG
#    		COMMANDS.append(CMD)

#logging.info("Start testing...")
#pool = Pool(16)  # four concurrent commands at a time
#for i, returncode in enumerate(pool.imap(partial(call, shell=True), COMMANDS)):
#	print('comando: ', COMMANDS[i])
#	logging.info('comando: ', COMMANDS[i])
#	if returncode != 0:
#		logging.info("%d command failed: %d" % (i, returncode))

for T in [10]:
	for dt in Gs:
		TEST = join(DIRTEST, dt)

		LOG = join('./tmp/{}.test.{}.log'.format(dt, T))
		CMD = './bin/mssRank ' + TEST + ' {} > '.format(T) + LOG
		COMMANDS.append(CMD)

logging.info("Start testing...")
pool = Pool(1)  # four concurrent commands at a time
for i, returncode in enumerate(pool.imap(partial(call, shell=True), COMMANDS)):
	print('comando: ', COMMANDS[i])
	logging.info('comando: ', COMMANDS[i])
	if returncode != 0:
		logging.info("%d command failed: %d" % (i, returncode))

