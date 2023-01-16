import os

base_folder = os.getcwd()

try:
    os.mkdir(os.path.expanduser('~/opt/anaconda3/envs'))
except:
    pass

os.chdir(os.path.expanduser('~/opt/anaconda3/envs'))
os.system('conda env create --file {}/pyrap_environment_alternative.yml --prefix pyRAP_trial'.format(base_folder))
