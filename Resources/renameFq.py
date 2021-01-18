import os
import glob
from glob import iglob
import random


os.chdir("/Users/rmvpaeme/test/BGZ_seq")

rootdir_glob = '/Users/rmvpaeme/test/BGZ_seq/**/*' # Note the added asterisks
# This will return absolute paths
file_list = [f for f in iglob('**/*', recursive=True) if os.path.isfile(f)]

for f in file_list:
    if f.startswith("RNA"):
        randomend = str(random.randrange(100000000, 200000000))
        foldername = f.split("/")[0]
        filename = f.split("/")[1]
        if not filename.startswith("RNA"):
            endname = '_'.join(filename.split("_")[1:])
            newname = foldername + "_" + endname
            newfullname = foldername + "/" + newname
            os.rename(f, newfullname)
