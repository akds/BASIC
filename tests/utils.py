import subprocess
import os


def find_bowtie2():
    """ Looks for bowtie2 in PATH if specific environmental variable does
        not exist
    """

    try:
        return os.environ['bowtie2']
    except KeyError:
        stdout, stderr = subprocess.Popen(["which", "bowtie2"],
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE).communicate()
        if stdout:
            return stdout.decode('ascii').strip()
        else:
            print('Bowtie2 not found. Please add it to $PATH or export the'
                  ' environmental variable "bowtie2" specifying its '
                  'directory.')
            exit(1)
