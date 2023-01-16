import os
from optparse import OptionParser
import sys
import glob
import shutil

parser = OptionParser('Usage:  -i input_folder -c chromsize_file -o output_folder')
parser.add_option("-i", "--input", dest="input",
                  help="Path to input folder", metavar="FOLDER")
parser.add_option("-c", "--chrom", dest="chrom",
                  help="Path to chrom_size file", metavar="FILE")
parser.add_option("-o", "--output", dest="output",
                  help="Path to output_folder", metavar="FOLDER")
(options, args) = parser.parse_args()

class ParserError(Exception):
    pass

if len(sys.argv) < 7:
    raise ParserError('''Not enough variables passed.
    'Usage: wigToBigWigBatch.py -i FOLDER -c FILE -o FOLDER'
    For more info type: wigToBigWigBatch.py --help''')

if options.input == options.output:
    raise ParserError('''Input folder can't be the same as output folder''')


if os.path.isdir(options.output) == True:
    os.system('rm -R {}'.format(options.output))

for file in glob.glob(options.input + '/*.gz'):
    os.system('gunzip {}'.format(file))


os.system('cp -R {} {}'.format(options.input, options.output))

### remove first line of wig file if it starts with 'track'
### this is necessary to run the bigwig converter wigToBigWig
for file in glob.glob(options.output + '/*.wig'):
    with open(file, 'r') as f_in:
        data = f_in.read().splitlines(True)
        f_in.close()
        n = 0
        if data[0].find('track') == 0:
            n += 1
        with open(file, 'w') as f_out:
            f_out.writelines(data[n:])
            f_out.close()


for file in glob.glob(options.output + '/*.wig'):
    os.system('wigToBigWig {} {} {}/{}.bw'.format(file, options.chrom, options.output, os.path.basename(file)))


for file in glob.glob(options.output + '/*.wig'):
    os.system('rm {}'.format(file))
