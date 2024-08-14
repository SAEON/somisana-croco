
import os
import subprocess

def conda_install(environment):
    subprocess.run(["conda", "env", "create", "-f", str(environment)],check=True)
