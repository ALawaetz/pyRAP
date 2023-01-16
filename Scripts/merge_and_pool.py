import os
import pandas as pd
import glob

__author__		= "Andreas Lawaetz"
__copyright__	= "Copyright 2022"
__version__		= "0.0.1"
__credits__		= ["Andreas Lawaetz"]
__maintainer__	= "Andreas Lawaetz"
__email__		= "acl58@bath.ac.uk"
__status__		= "Production"

##################################################################################
#
#	options_unique.py
#
#
#	Copyright (c) Andreas Lawaetz 2022
#
#	Permission is hereby granted, free of charge, to any person obtaining a copy
#	of this software and associated documentation files (the "Software"), to deal
#	in the Software without restriction, including without limitation the rights
#	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#	copies of the Software, and to permit persons to whom the Software is
#	furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#	THE SOFTWARE.
#
##################################################################################

# # ############### Merge forward and reverse annotations #######################

def merge_and_pool(gff_folder, output_folder, list_of_conditions):
    tracklist = []
    for file in sorted(glob.glob(gff_folder + '*.gff3')):
        tracklist.append(os.path.basename(file))

    for i in range(0, len(tracklist), 2):
        df = pd.read_csv(gff_folder + tracklist[i], sep = '\t', header = 0)
        df = pd.concat([df, pd.read_csv(gff_folder + tracklist[i + 1], sep = '\t', header = 0)])
        df.to_csv(os.path.join(output_folder, 'FwdandRev_' + tracklist[i]), sep = '\t', index = False)
        os.remove(gff_folder + tracklist[i])
        os.remove(gff_folder + tracklist[i + 1])
        del df


    # # ### In attributes we write in what condition any transcript is expressed using
    # ### the list of conditions found at the beginning of the pipeline.

    i = 0
    for file in sorted(glob.glob(output_folder + 'FwdandRev_*')):
        df = pd.read_csv(file, sep = '\t', header = 0)
        attributes = [list_of_conditions[i]] * len(df)
        df['attributes'] = attributes
        df.to_csv(os.path.join(output_folder, list_of_conditions[i] + '.gff3'), sep = '\t', index = False)
        os.remove(file)
        i += 1

    #
    # ################### Pool all annotations ##############################
    #
    pool = pd.DataFrame()

    for file in sorted(glob.glob(output_folder + '*.gff3')):
        df = pd.read_csv(file, sep = '\t', header = 0)
        pool = pd.concat([pool, df])

    pool = pool.sort_values(['start', 'end'], axis = 0, ascending = True)
    pool.to_csv(os.path.join(output_folder, 'All_conditions.gff3'), sep = '\t', index = False)

def pool(folder):
    pool = pd.DataFrame()

    for file in sorted(glob.glob(folder + '*.gff3')):
        df = pd.read_csv(file, sep = '\t', header = 0)
        pool = pd.concat([pool, df])

    pool = pool.sort_values(['start', 'end'], axis = 0, ascending = True)
    pool.to_csv(os.path.join(folder, 'All_conditions.gff3'), sep = '\t', index = False)
