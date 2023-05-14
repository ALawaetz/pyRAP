import os, shutil
## for repeated run through of the pipeline without the need to manually
# delete folders, the following 6 lines are in place.

try:
    os.chdir('Pipeline')
    if '.DS_Store' in os.listdir():
        os.remove('.DS_Store')
    for folder in os.listdir():
        shutil.rmtree(folder, ignore_errors=False, onerror=None)
    os.chdir('..')
except:
    pass
