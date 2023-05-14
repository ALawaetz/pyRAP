import sys
import pandas as pd
import os
from optparse import OptionParser
import sys
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

parser = OptionParser('''Usage: Read_depth.py --folder path_to_wig_folder --ex condition_name --size path_to_chrom_size_file
Calulate read depth in reads pr million.
Take the sum of all reads in the 4 wig files constituting a condition
Divide by 4
Divide by genome size''')
parser.add_option("--folder", dest="wig_folder",
                  help="folder", metavar="folder")
parser.add_option("--ex", dest="condition",
                  help="condition_name", metavar="string")
parser.add_option("--size", dest="chrom",
                  help="path to chrom size file", metavar="file")                 
(options, args) = parser.parse_args()


class ParserError(Exception):
    pass
if len(sys.argv) < 3:
    raise ParserError('''Not passed enough variables
    Usage: Usage: Read_depth.py --folder path_to_wig_folder --ex condition_name --size path_to_chrom_size_file
    For more info type: Read_depth.py --help''')

thesum = []
for file in glob.glob(options.wig_folder + '/*{}*'.format(options.condition)):
    with open(file, 'r') as f:
        text = f.read().splitlines(True)
        f.close()
        n = 0
        for item in text:
            alist = item.split('\t')
            try:
                float(alist[0])
                break
            except ValueError:
                n += 1

    df = pd.read_csv(file, sep = '\t', names = ('location', 'value'), skiprows = n)
    thesum.append(sum(df['value'].tolist()))

sum_of_four = sum(thesum)

chrom = pd.read_csv(options.chrom, sep = '\t', names = ['chromosome', 'size'])
genome_size = chrom['size'][0]

mean_read_value = sum_of_four / 4 / genome_size

print('Mean read value for condition: {} is {} reads pr. base'.format(options.condition, mean_read_value))
