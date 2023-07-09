import os

base_folder = os.getcwd()

try:
    os.mkdir(os.path.expanduser('~/anaconda3/envs'))
except:
    pass

os.chdir(os.path.expanduser('~/anaconda3/envs'))
os.system('conda env create --file {}/pyrap_environment.yml --prefix pyRAP'.format(base_folder))
