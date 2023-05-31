import psutil
import subprocess
import time
import os

class ProcessTimer:
	"""Class is a modification of https://stackoverflow.com/questions/13607391/subprocess-memory-usage-in-python"""
	def __init__(self,command,inp):
		self.command = command
		self.input = inp
		self.execution_state = False

	def execute(self):
		self.max_vms_memory = 0
		self.max_rss_memory = 0

		self.p = subprocess.Popen(args=self.command,shell=False, stdin=subprocess.PIPE)
		self.p.stdin.write(bytes(self.input[0] + "\n", 'utf-8'))
		self.p.stdin.flush()
		time.sleep(0.1)
		self.p.stdin.write(bytes(self.input[1] + "\n", 'utf-8'))
		self.p.stdin.flush()
		self.t1 = None
		self.t0 = time.time()
		self.execution_state = True

	def poll(self):
		if not self.check_execution_state():
			return False
		
		self.t1 = time.time()

		try:
			pp = psutil.Process(self.p.pid)

			descendants = [pp]

			rss_memory = 0
			vms_memory = 0

			for descendant in descendants:
				try:
					mem_info = descendant.memory_info()
					rss_memory += mem_info[0]
					vms_memory += mem_info[1]
				except psutil.NoSuchProcess:
					pass
			self.max_vms_memory = max(self.max_vms_memory,vms_memory)
			self.max_rss_memory = max(self.max_rss_memory,rss_memory)

		except psutil.NoSuchProcess:
			return self.check_execution_state()

		return self.check_execution_state()


	def is_running(self):
		return psutil.pid_exists(self.p.pid) and self.p.poll() == None
	def check_execution_state(self):
		if not self.execution_state:
			return False
		if self.is_running():
			return True
		self.executation_state = False
		self.t1 = time.time()
		return False

	def close(self,kill=False):
		try:
			pp = psutil.Process(self.p.pid)
			if kill:
				pp.kill()
			else:
				pp.terminate()
		except psutil.NoSuchProcess:
			pass

EXAMPLES_DIR = 'data/'

def main():
		for dir in os.listdir(EXAMPLES_DIR):
				if not dir.startswith('J31'):
						continue
				
				filename = EXAMPLES_DIR + dir
				ptimer = ProcessTimer(["build/src/GVFinder"], [filename, "out.txt"])
				try:
						ptimer.execute()
						while ptimer.poll():
								time.sleep(.1)
				finally:
						ptimer.close()

				print(filename)		
				print ('return code:',ptimer.p.returncode)
				print ('time:',ptimer.t1 - ptimer.t0, "seconds")
				print ('max virtual memory usage:',ptimer.max_vms_memory/1024, "KiB")
				print ('max physical memory usage:',ptimer.max_rss_memory/1024, "KiB")
				
if __name__ == '__main__':
	main()