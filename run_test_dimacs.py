from os import listdir
from os.path import isfile, join, isdir
from subprocess import check_output

from functools import partial
from multiprocessing.dummy import Pool
from subprocess import call

import logging

Gs = ['brock200_1.col','brock200_2.col','brock200_3.col','brock200_4.col',
   'brock400_1.col','brock400_2.col','brock400_3.col','brock400_4.col',
   'C125.9.col', 'C250.9.col', 'C500.9.col','keller4.col', 'hamming6-4.col', 
   'p_hat300-2.col', 'p_hat300-3.col', 'MANN_09.col',
   'sanr200_0.7.col', 'sanr200_0.9.col', 'sanr400_0.7.col']
#'p_hat500-3.col',

# Gx TEST
DIRTEST = './data/dimacs/'
COMMANDS = []

#for T in [11, 12, 13, 14, 5]:
	#for dt in Gs:
    		#TEST = join(DIRTEST, dt)

    		#LOG = join('./tmp/{}.test.{}.log'.format(dt, T))
    		#CMD = './bin/mssRank ' + TEST + ' {} 1 > '.format(T) + LOG
    		#COMMANDS.append(CMD)

#logging.info("Start testing...")
#pool = Pool(16)  # four concurrent commands at a time
#for i, returncode in enumerate(pool.imap(partial(call, shell=True), COMMANDS)):
	#print('comando: ', COMMANDS[i])
	#logging.info('comando: ', COMMANDS[i])
	#if returncode != 0:
		#logging.info("%d command failed: %d" % (i, returncode))

for T in [10]:
	for dt in Gs:
		TEST = join(DIRTEST, dt)

		LOG = join('./tmp/{}.test.{}.log'.format(dt, T))
		CMD = './bin/mssRank ' + TEST + ' {} 1 > '.format(T) + LOG
		COMMANDS.append(CMD)

logging.info("Start testing...")
pool = Pool(1)  # four concurrent commands at a time
for i, returncode in enumerate(pool.imap(partial(call, shell=True), COMMANDS)):
	print('comando: ', COMMANDS[i])
	logging.info('comando: ', COMMANDS[i])
	if returncode != 0:
		logging.info("%d command failed: %d" % (i, returncode))

