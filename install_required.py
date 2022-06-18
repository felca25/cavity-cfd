import sys
import subprocess


# implementing pip as a subprocess:
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numpy'])
# process output with an API in the subprocess module:
reqs = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
installed_packages = [r.decode().split('==')[0] for r in reqs.split()]
print(installed_packages)