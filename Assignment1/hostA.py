import os
import subprocess

dir = os.getcwd()
#print(dir)
dir_list =  os.listdir(dir)
#print(dir_list)

# deleting already exist hostfile
if 'hostsfile' in dir_list:
	os.system('rm hostsfile')

# checking host in cluster
for host in range(1,41):
	try:
    		response = subprocess.check_output(
        		['ping', '-c1', '-w1', 'csews'+str(host)],
        		stderr=subprocess.STDOUT,  # get all output
        		universal_newlines=True  # return string not bytes
    			)
		with open('hostsfile','a') as file:
			file.write('csews'+str(host)+":8\n")

	except subprocess.CalledProcessError:
    		response = None

print("Hostfile Created!!!")
