import install_required as inst_req

import sys
import subprocess
import conda.cli.python_api as Conda


# implementing pip as a subprocess:
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numpy'])
# process output with an API in the subprocess module:
reqs = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
installed_packages = [r.decode().split('==')[0] for r in reqs.split()]
print(installed_packages)

def conda_install():
    ['os', 'json', 'time', 'numpy', 'matplotlib', 'numba',]

if __name__ == "__main__":
    