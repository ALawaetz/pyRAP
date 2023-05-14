__author__		= "Andreas Lawaetz"
__copyright__	= "Copyright 2022"
__version__		= "0.0.1"
__credits__		= ["Andreas Lawaetz"]
__maintainer__	= "Andreas Lawaetz"
__email__		= "acl58@bath.ac.uk"
__status__		= "Production"

##################################################################################
#
#	GTF_to_GFF.py
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

import sys
import pandas as pd
import os
from optparse import OptionParser




parser = OptionParser('Usage: GTF_to_GFF.py -f inputfile -o outputfile')
parser.add_option("-f", "--file", dest="input_file",
                  help="input file", metavar="FILE")
parser.add_option("-o", "--output", dest="output_file",
                  help="output file", metavar="FILE")
(options, args) = parser.parse_args()


class ParserError(Exception):
    pass
if len(sys.argv) < 3:
    raise ParserError('''Not passed enough variables
    Usage: GTF_to_GFF.py input output
    For more info type: GTF_to_GFF.py --help''')


with open(options.input_file, 'r') as f:
    alist = f.readline().split('\t')
    try:
        header = int(float(alist[3]))
    except ValueError:
        header = 'yes'
    if type(header) == int:
        df = pd.read_csv(options.input_file, sep = '\t', names = ['seqID','source','feature','start','end','score','strand','phase', 'attributes'])
    else:
        df = pd.read_csv(options.input_file, sep = '\t', header = 0)

attributes = df['attributes'].tolist()

new_attributes = []
for item in attributes:
    new_attributes.append(item.replace(' ', '='))

df['attributes'] = new_attributes

thepath = options.output_file.replace(os.path.basename(options.output_file), '')

df.to_csv(os.path.join(thepath, os.path.basename(options.output_file) + '.gff3'), sep = '\t', header = None, index = False)
