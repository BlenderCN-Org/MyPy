#!/usr/bin/ipython console

#class CreateConnection(object):

import paramiko

client = paramiko.SSHClient()
client.load_system_host_keys()
client.set_missing_host_key_policy(paramiko.WarningPolicy())

client.connect('128.178.134.191',username='riccardo',password='06111983')
#stdin,stdout,stderr = client.exec_command('ls')


#client.close()
