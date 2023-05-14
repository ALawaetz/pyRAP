import sys
import pandas as pd
import os
from optparse import OptionParser
import sys


__author__		= "Andreas Lawaetz"
__copyright__	= "Copyright 2022"
__version__		= "0.0.1"
__credits__		= ["Andreas Lawaetz"]
__maintainer__	= "Andreas Lawaetz"
__email__		= "acl58@bath.ac.uk"
__status__		= "Production"

##################################################################################
#
#	GFF_to_GTF.py
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

parser = OptionParser('''Usage: Overlap.py --f1 file1 --file2 file2
Finds annnotations in file1 not overlapping annotations in file2 and write to output file
output file is named by contract the file names of file1 and file2 and the file is saved
in the same folder as file1''')
parser.add_option("--f1", dest="file1",
                  help="file1", metavar="FILE")
parser.add_option("--f2", dest="file2",
                  help="file2", metavar="FILE")                 
(options, args) = parser.parse_args()


class ParserError(Exception):
    pass
if len(sys.argv) < 3:
    raise ParserError('''Not passed enough variables
    Usage: Usage: Overlap.py --f1 file1 --file2 file2
    For more info type: Overlap.py --help''')


with open(options.file1, 'r') as f:
    alist = f.readline().split('\t')
    try:
        header = int(float(alist[3]))
    except ValueError:
        header = 'yes'
    if type(header) == int:
        df1 = pd.read_csv(options.file1, sep = '\t', names = ['seqID','source','feature','start','end','score','strand','phase', 'attributes'])
    else:
        df1 = pd.read_csv(options.file1, sep = '\t', header = 0)

with open(options.file2, 'r') as f:
    alist = f.readline().split('\t')
    try:
        header = int(float(alist[3]))
    except ValueError:
        header = 'yes'
    if type(header) == int:
        df2 = pd.read_csv(options.file2, sep = '\t', names = ['seqID','source','feature','start','end','score','strand','phase', 'attributes'])
    else:
        df2 = pd.read_csv(options.file2, sep = '\t', header = 0)


non_overlap = pd.DataFrame()
for index, row in df1.iterrows():
    R_strand = row['strand']
    R_start = row['start']
    R_end = row['end']
    ### If there is an overlap then one of these conditions will provide a dataframe in df2 with a length above 0
    conditions = (R_start <= df2['start']) & (R_end > df2['start']) | (R_start >= df2['start']) & (R_start < df2['end']) | (R_start == df2['start']) & (R_end == df2['end'])
    if len(df2[conditions][df2[conditions]['strand'] == R_strand]) == 0:
        non_overlap = pd.concat([non_overlap, df1[df1.index == index]])


non_overlap.to_csv(options.file1 + '_' + os.path.basename(options.file2) + '_no_overlap.gff3', sep = '\t', index = False, header = None)