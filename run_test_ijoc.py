from os import listdir
from os.path import isfile, join, isdir
from subprocess import check_output

from functools import partial
from multiprocessing.dummy import Pool
from subprocess import call

import logging

DIRTEST = 'C:\\Users\\Gualandi\\Desktop\\data-rudy\\'
COMMANDS = []

Gs = [
    'r-50-0.1-1-14', 'r-50-0.1-2-15', 'r-50-0.1-3-16', 'r-50-0.1-4-17',
    'r-50-0.1-5-18', 'r-50-0.1-6-19', 'r-50-0.1-7-20', 'r-50-0.1-8-21',
    'r-50-0.1-9-22', 'r-50-0.1-10-23', 'r-50-0.2-1-24', 'r-50-0.2-2-25',
    'r-50-0.2-3-26', 'r-50-0.2-4-27', 'r-50-0.2-5-28', 'r-50-0.2-6-29',
    'r-50-0.2-7-30', 'r-50-0.2-8-31', 'r-50-0.2-9-32', 'r-50-0.2-10-33',
    'r-50-0.3-1-34', 'r-50-0.3-2-35', 'r-50-0.3-3-36', 'r-50-0.3-4-37',
    'r-50-0.3-5-38', 'r-50-0.3-6-39', 'r-50-0.3-7-40', 'r-50-0.3-8-41',
    'r-50-0.3-9-42', 'r-50-0.3-10-43', 'r-50-0.4-1-44', 'r-50-0.4-2-45',
    'r-50-0.4-3-46', 'r-50-0.4-4-47', 'r-50-0.4-5-48', 'r-50-0.4-6-49',
    'r-50-0.4-7-50', 'r-50-0.4-8-51', 'r-50-0.4-9-52', 'r-50-0.4-10-53',
    'r-50-0.5-1-54', 'r-50-0.5-2-55', 'r-50-0.5-3-56', 'r-50-0.5-4-57',
    'r-50-0.5-5-58', 'r-50-0.5-6-59', 'r-50-0.5-7-60', 'r-50-0.5-8-61',
    'r-50-0.5-9-62', 'r-50-0.5-10-63', 'r-50-0.6-1-64', 'r-50-0.6-2-65',
    'r-50-0.6-3-66', 'r-50-0.6-4-67', 'r-50-0.6-5-68', 'r-50-0.6-6-69',
    'r-50-0.6-7-70', 'r-50-0.6-8-71', 'r-50-0.6-9-72', 'r-50-0.6-10-73',
    'r-50-0.7-1-74', 'r-50-0.7-2-75', 'r-50-0.7-3-76', 'r-50-0.7-4-77',
    'r-50-0.7-5-78', 'r-50-0.7-6-79', 'r-50-0.7-7-80', 'r-50-0.7-8-81',
    'r-50-0.7-9-82', 'r-50-0.7-10-83', 'r-50-0.8-1-84', 'r-50-0.8-2-85',
    'r-50-0.8-3-86', 'r-50-0.8-4-87', 'r-50-0.8-5-88', 'r-50-0.8-6-89',
    'r-50-0.8-7-90', 'r-50-0.8-8-91', 'r-50-0.8-9-92', 'r-50-0.8-10-93',
    'r-50-0.9-1-94', 'r-50-0.9-2-95', 'r-50-0.9-3-96', 'r-50-0.9-4-97',
    'r-50-0.9-5-98', 'r-50-0.9-6-99', 'r-50-0.9-7-100', 'r-50-0.9-8-101',
    'r-50-0.9-9-102', 'r-50-0.9-10-103', 'r-75-0.1-1-104', 'r-75-0.1-2-105',
    'r-75-0.1-3-106', 'r-75-0.1-4-107', 'r-75-0.1-5-108', 'r-75-0.1-6-109',
    'r-75-0.1-7-110', 'r-75-0.1-8-111', 'r-75-0.1-9-112', 'r-75-0.1-10-113',
    'r-75-0.2-1-114', 'r-75-0.2-2-115', 'r-75-0.2-3-116', 'r-75-0.2-4-117',
    'r-75-0.2-5-118', 'r-75-0.2-6-119', 'r-75-0.2-7-120', 'r-75-0.2-8-121',
    'r-75-0.2-9-122', 'r-75-0.2-10-123', 'r-75-0.3-1-124', 'r-75-0.3-2-125',
    'r-75-0.3-3-126', 'r-75-0.3-4-127', 'r-75-0.3-5-128', 'r-75-0.3-6-129',
    'r-75-0.3-7-130', 'r-75-0.3-8-131', 'r-75-0.3-9-132', 'r-75-0.3-10-133',
    'r-75-0.4-1-134', 'r-75-0.4-2-135', 'r-75-0.4-3-136', 'r-75-0.4-4-137',
    'r-75-0.4-5-138', 'r-75-0.4-6-139', 'r-75-0.4-7-140', 'r-75-0.4-8-141',
    'r-75-0.4-9-142', 'r-75-0.4-10-143', 'r-75-0.5-1-144', 'r-75-0.5-2-145',
    'r-75-0.5-3-146', 'r-75-0.5-4-147', 'r-75-0.5-5-148', 'r-75-0.5-6-149',
    'r-75-0.5-7-150', 'r-75-0.5-8-151', 'r-75-0.5-9-152', 'r-75-0.5-10-153',
    'r-75-0.6-1-154', 'r-75-0.6-2-155', 'r-75-0.6-3-156', 'r-75-0.6-4-157',
    'r-75-0.6-5-158', 'r-75-0.6-6-159', 'r-75-0.6-7-160', 'r-75-0.6-8-161',
    'r-75-0.6-9-162', 'r-75-0.6-10-163', 'r-75-0.7-1-164', 'r-75-0.7-2-165',
    'r-75-0.7-3-166', 'r-75-0.7-4-167', 'r-75-0.7-5-168', 'r-75-0.7-6-169',
    'r-75-0.7-7-170', 'r-75-0.7-8-171', 'r-75-0.7-9-172', 'r-75-0.7-10-173',
    'r-75-0.8-1-174', 'r-75-0.8-2-175', 'r-75-0.8-3-176', 'r-75-0.8-4-177',
    'r-75-0.8-5-178', 'r-75-0.8-6-179', 'r-75-0.8-7-180', 'r-75-0.8-8-181',
    'r-75-0.8-9-182', 'r-75-0.8-10-183', 'r-75-0.9-1-184', 'r-75-0.9-2-185',
    'r-75-0.9-3-186', 'r-75-0.9-4-187', 'r-75-0.9-5-188', 'r-75-0.9-6-189',
    'r-75-0.9-7-190', 'r-75-0.9-8-191', 'r-75-0.9-9-192', 'r-75-0.9-10-193',
    'r-100-0.1-1-194', 'r-100-0.1-2-195', 'r-100-0.1-3-196', 'r-100-0.1-4-197',
    'r-100-0.1-5-198', 'r-100-0.1-6-199', 'r-100-0.1-7-200', 'r-100-0.1-8-201',
    'r-100-0.1-9-202', 'r-100-0.1-10-203', 'r-100-0.2-1-204', 'r-100-0.2-2-205',
    'r-100-0.2-3-206', 'r-100-0.2-4-207', 'r-100-0.2-5-208', 'r-100-0.2-6-209',
    'r-100-0.2-7-210', 'r-100-0.2-8-211', 'r-100-0.2-9-212', 'r-100-0.2-10-213',
    'r-100-0.3-1-214', 'r-100-0.3-2-215', 'r-100-0.3-3-216', 'r-100-0.3-4-217',
    'r-100-0.3-5-218', 'r-100-0.3-6-219', 'r-100-0.3-7-220', 'r-100-0.3-8-221',
    'r-100-0.3-9-222', 'r-100-0.3-10-223', 'r-100-0.4-1-224', 'r-100-0.4-2-225',
    'r-100-0.4-3-226', 'r-100-0.4-4-227', 'r-100-0.4-5-228', 'r-100-0.4-6-229',
    'r-100-0.4-7-230', 'r-100-0.4-8-231', 'r-100-0.4-9-232', 'r-100-0.4-10-233',
    'r-100-0.5-1-234', 'r-100-0.5-2-235', 'r-100-0.5-3-236', 'r-100-0.5-4-237',
    'r-100-0.5-5-238', 'r-100-0.5-6-239', 'r-100-0.5-7-240', 'r-100-0.5-8-241',
    'r-100-0.5-9-242', 'r-100-0.5-10-243', 'r-100-0.6-1-244', 'r-100-0.6-2-245',
    'r-100-0.6-3-246', 'r-100-0.6-4-247', 'r-100-0.6-5-248', 'r-100-0.6-6-249',
    'r-100-0.6-7-250', 'r-100-0.6-8-251', 'r-100-0.6-9-252', 'r-100-0.6-10-253',
    'r-100-0.7-1-254', 'r-100-0.7-2-255', 'r-100-0.7-3-256', 'r-100-0.7-4-257',
    'r-100-0.7-5-258', 'r-100-0.7-6-259', 'r-100-0.7-7-260', 'r-100-0.7-8-261',
    'r-100-0.7-9-262', 'r-100-0.7-10-263', 'r-100-0.8-1-264', 'r-100-0.8-2-265',
    'r-100-0.8-3-266', 'r-100-0.8-4-267', 'r-100-0.8-5-268', 'r-100-0.8-6-269',
    'r-100-0.8-7-270', 'r-100-0.8-8-271', 'r-100-0.8-9-272', 'r-100-0.8-10-273',
    'r-100-0.9-1-274', 'r-100-0.9-2-275', 'r-100-0.9-3-276', 'r-100-0.9-4-277',
    'r-100-0.9-5-278', 'r-100-0.9-6-279', 'r-100-0.9-7-280', 'r-100-0.9-8-281',
    'r-100-0.9-9-282', 'r-100-0.9-10-283' ] 

DIRTEST = './data/rnd/'
COMMANDS = []

for T in [10]:
	for dt in Gs:
		TEST = join(DIRTEST, dt)

		LOG = join('./tmp/{}.test.{}.log'.format(dt, T))
		CMD = './bin/mssRank ' + TEST + ' {} > '.format(T) + LOG
		COMMANDS.append(CMD)

logging.info("Start testing...")
print(COMMANDS)
pool = Pool(10)  # four concurrent commands at a time
for i, returncode in enumerate(pool.imap(partial(call, shell=True), COMMANDS)):
    print('comando: ', COMMANDS[i])
    if returncode != 0:
        print("%d command failed: %d" % (i, returncode))
