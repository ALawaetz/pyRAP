__author__		= "Andreas Lawaetz"
__copyright__	= "Copyright 2022"
__version__		= "0.0.1"
__credits__		= ["Andreas Lawaetz"]
__maintainer__	= "Andreas Lawaetz"
__email__		= "acl58@bath.ac.uk"
__status__		= "Production"

##################################################################################
#
#	1_Rend_seq_pipeline_22August2022.py
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

import pandas as pd
import numpy as np
import os
import glob
from multiprocessing import Process, Pool
import subprocess
from optparse import OptionParser
import sys
import shutil
from Scripts import folders, merge_and_pool, options_unique, pool_within, post_processing, statistical_analysis, Draw_GFF3_supply_w_capterm
from Scripts import repeated_runs

parser = OptionParser('Usage: 1_Rend_seq_pipeline_22August.py [options] -f inputfolder')
parser.add_option("-f", "--input_folder", dest="wig_folder",
                  help="provide the path to the folder containing your wig files", metavar="FOLDER")
parser.add_option("-n", "--names_file", dest="names_file",
                  help="Location of txt file containing names of wig files (conditions). Default is file names in wig folder. List of names has to be provided in the same order as the order of files in the wig folder when it is alphabetically sorted.", metavar="FILE")
parser.add_option("-c", "--chrom", dest="chrom",
                  help="path to tab-seperated file containing chromosome name and size.", metavar="file")
parser.add_option("--genome_start", dest="genome_start", type="int",
                  help="Provide start coordinate if only a subset of genome wants analysed", metavar="integer", default=0)
parser.add_option("--genome_end", dest="genome_end", type="int",
                  help="Provide end coordinate if only a subset of genome wants analysed", metavar="integer")
parser.add_option("-g", "--guide", dest="operon_guide",
                  help="provide name of file used to define operons. Name has to match file name in wig folder or name in list of names if provided (-n flag).", metavar="wig_file_name")
parser.add_option("--operon_cutoff", dest="operon_cutoff", type='float',
                  help="Read density to define operons. Choose value between 0 and 1.", metavar="float", default=0.8)
parser.add_option("-o", "--operons", dest="operons",
                  help="path to operon file in GFF format. File shall contain no header.", metavar="file")
parser.add_option("--alternative_operons", dest="alternative_operons",
                  help="path to alternative operon file in GFF or GTF file format", metavar="wig_file_name")
parser.add_option("--z_start", dest="z_start", type="float",
                  help="z_score to call start peaks", metavar="float", default=7)
parser.add_option("--z_stop", dest="z_stop", type="float",
                  help="z_score to call stop peaks", metavar="float", default=6)
parser.add_option("--z_small", dest="z_small", type="float",
                  help="z_score to call peaks in small windows", metavar="float", default=5)
parser.add_option("--z_ii", dest="z_ii", type="float",
                  help="z_score to reevalute unique peaks.", metavar="float", default=7)
parser.add_option("--unique", dest="unique", type='int',
                  help="This flag can be choosen to reevaluate peaks only found in N number of samples. E.g. --unique 2 will reevaluate peaks only found in 2 or less samples. Z score to reevaluate peaks is determined with z_ii flag.", metavar="integer")
parser.add_option("--fdr", dest="fdr", type="float",
                  help="FDR z_score to reevaluate low-confidence peaks. Default is 2. Choose 0 to surpass", metavar="integer", default=2)
parser.add_option("--old_annotation", dest="old_annotation",
                  help="Path to old annotation file. Used to transfer names and features when overlap and for statistical analysis", metavar="FILE")
parser.add_option("--sRNA_features", dest="sRNA_features",
                  help="Choose this option to define what features in old annotation file is sRNAs. Default = S_feature, non-coding RNA, misc_RNA, ncRNA, antisense_RNA", metavar='string', type='str', default = 'S_feature, non-coding RNA, misc_RNA, ncRNA, antisense_RNA')
parser.add_option("--source_exclude", dest="source_exclude",
                  help="If old annotation is based on different sources this flag can be used to exclude one of them in the statistical analysis. For example setting the flag to 'prediction' will exclude annotations where the source (column 2 in GFF3 file format) is prediction", metavar="STRING")
parser.add_option("--format", dest="format", choices=['GFF', 'GTF'],
                  help="Specify file format of old_annotation file (GFF or GTF)", default = 'GFF')
parser.add_option("--pool", dest="pool", type = 'int',
                  help="Pool peaks that are within N bp of each other", metavar="integer", default=5)
parser.add_option("--output_folder", dest="output_folder",
                  help="Path to output folder", metavar="FOLDER")
parser.add_option("--output", dest="output", choices=['GFF3', 'GTF', 'All'],
                  help="Define outout (GFF3, GTF, or All)", default='GFF3')
parser.add_option("--statistics", dest="statistics", action='store_true',
                  help="Choose this option if you want to add a statistical analysis")
parser.add_option("--stat_cond", dest="stat_cond",
                  help="Choose condition to do statistics on. Default is all conditions", default='All_conditions', metavar='string', type='str')
parser.add_option("--compare_starts", dest="cap_seq",
                  help="Path to start sites from alternative method. Has to be comma seperated csv file.", metavar='FILE')
parser.add_option("--compare_stops", dest="term_seq",
                  help="Path to stop sites from alternative method. Has to be comma seperated csv file.", metavar='FILE')
parser.add_option("--search_distance", dest="search_distance",
                  help="Defines distance to verify with flags --compare_starts and --compare_stops. Default=5", metavar='integer', type='int', default=5)
parser.add_option("--operon_starts", dest="operon_starts",
                  help="Genomic positions of operon start sites", metavar='FILE')
parser.add_option("--operon_stops", dest="operon_stops",
                  help="Genomic positions of operon stop sites", metavar='FILE')
parser.add_option("--search_distance_0", dest="search_distance_0",
                  help="Defines distance to verify with flags --operon_starts and --operon_stops. Default=30", metavar='integer', type='int', default=30)
parser.add_option("--compare_condition", dest="compare_condition",
                  help="Condition to compare to Cap-seq and Term-seq", metavar='string', default='All_conditions')
parser.add_option("--feature_distance", dest="feature_distance",
                  help="Provide a txt file to create histogram plots showing the nucleotide distance between an old annotation file and the new Rend-seq annotation file. File should be tab seperated with 4 columns and a header. Column names are feature, dist_id, sRNA_type and title.", metavar='FILE', default=None)
parser.add_option("--compare_peaks_across_samples", dest="compare_peaks_across_samples",
                  help="Provide a txt file with the information of what conditions to compare. File should contain two columns and a header. The name of the columns are peaks_in and peaks_not_in. Names should match names given with the -n flag.", metavar='FILE', default=None)
parser.add_option("--z_min_compare", dest="z_min_compare",
                  help="Minimum z value to to pass peaks for comparison. Use with --compare_peaks_across_samples flag", metavar='INTEGER', type='int', default=None)
parser.add_option("--z_max_compare", dest="z_max_compare",
                  help="Maximum z value to to pass non-peaks for comparison. Use with --compare_peaks_across_samples flag", metavar='INTEGER', type='int', default=None)
parser.add_option("--density_compare", dest="density_compare",
                  help="Read density cutoff for differentially expressed genes determination. Use with --compare_peaks_across_samples flag", metavar='FLOAT', type='float', default=0.5)
parser.add_option("-p", "--processors", dest="processors",
                  help="Number of processors to use. Default is 8", metavar='INT', type='int', default=8)


(options, args) = parser.parse_args()


class ParserError(Exception):
    pass
if options.wig_folder == None:
    raise ParserError('''No inputfolder provided.
    Usage: 1_Rend_seq_pipeline_22August.py -f FOLDER -c FILE --output_folder FOLDER
    For more info type: 1_Rend_seq_pipeline_22August.py --help''')

if options.chrom == None:
    raise ParserError('''No chrom file provided.
    Usage: 1_Rend_seq_pipeline_22August.py -f FOLDER -c FILE --output_folder FOLDER
    For more info type: 1_Rend_seq_pipeline_22August.py --help''')

if options.output_folder == None:
    raise ParserError('''No output_folder specified.
    Usage: 1_Rend_seq_pipeline_22August.py -f FOLDER -c FILE --output_folder FOLDER
    For more info type: 1_Rend_seq_pipeline_22August.py --help''')

if options.operon_cutoff < 0 or options.operon_cutoff > 1:
    raise ParserError('''Invalid operon cutoff. Choose value between 0 and 1.''')


chrom = pd.read_csv(options.chrom, sep = '\t', names = ['chromosome', 'size'])
genome_name = chrom['chromosome'][0]
genome_size = chrom['size'][0]
if len(chrom) != 1:
    raise ParserError('Chrom file must contain only one chromsome')

if options.genome_end == None:
    options.genome_end = genome_size

if options.genome_start == None or options.genome_end == None:
    genome_size = options.genome_end - options.genome_start

### if wig folder contain gzipped files then unzip
for file in glob.glob(options.wig_folder + '/*.gz'):
    os.system('gunzip {}'.format(file))

### List of file names
if options.names_file == None:
    list_of_conditions = []
    n = 0
    for file in sorted(glob.glob(options.wig_folder + '/*.wig')):
        if n % 4 == 0:
            list_of_conditions.append(os.path.basename(file)[:-4])
        else:
            pass
        n += 1
else:
    with open(options.names_file, 'r') as f:
        astring = f.read()
        raw_list_of_conditions = astring.split('\n')
        list_of_conditions = []
        for item in raw_list_of_conditions:
            if len(item) == 0:
                pass
            else:
                list_of_conditions.append(item)


if options.operon_guide == None:
    options.operon_guide = list_of_conditions[0]


if options.old_annotation != None:
    old_annotation_name = os.path.basename(options.old_annotation)


### Make copy of provided annotation file
annotation_file_list = []
if options.old_annotation != None:
    annotation_file_list.append(options.old_annotation)

for item in annotation_file_list:
    try:
        shutil.copy(item, folders.base_folder + '/{}'.format(os.path.basename(item)))
    except shutil.SameFileError:
        pass 
        

if options.old_annotation != None:
    options.old_annotation = folders.base_folder + '/' + os.path.basename(options.old_annotation)


annotation_file_list = []
if options.old_annotation != None:
    annotation_file_list.append(options.old_annotation)


for item in annotation_file_list:
    if options.format == 'GTF':
        os.system('python3 {} -f {} -o {}'.format(folders.GTF_to_GFF_script, item, item))
        if os.path.basename(item) + '.gff3' in os.listdir(folders.base_folder):
            pass
        else:
            os.wait()

### check if old annotation file provided starts with lines of #
### Annotation files from genbank start with an annotation that spans the entire genome. We remove that.
# if item != None:
    with open(item, 'r') as f_in:
        text = f_in.read().splitlines(True)
        f_in.close()
        n = 0
        for t in text:
            if t.find('#') == 0:
                n += 1
            else:
                break
        alist = text[n].split('\t')
        try:
            float(alist[3])
        except:
            n += 1
            alist = text[n].split('\t')
        if int(alist[3]) == 1 and int(alist[4]) == genome_size:
            n += 1
        else:
            pass
        with open(item, 'w') as f_out:
            f_out.writelines(text[n:])
            f_out.close()
# else:
#     pass


### check if old annotation file provided has header or not. If not, make one. Save in base_folder
# if item != None:
    with open(item, 'r') as f:
        alist = f.readline().split('\t')
        f.close()
        try:
            float(alist[3])
            df = pd.read_csv(item, sep = '\t', names = ['seqID','source','feature','start','end','score','strand','phase', 'attributes'])
            df.to_csv(os.path.join(folders.base_folder, os.path.basename(item)), sep = '\t', index = False)
        except ValueError:
            df = pd.read_csv(item, sep = '\t', skiprows = 1, names = ['seqID','source','feature','start','end','score','strand','phase', 'attributes'])
            df.to_csv(os.path.join(folders.base_folder, os.path.basename(item)), sep = '\t', index = False)


if options.z_min_compare == None:
    options.z_min_compare = options.z_start + 1

if options.z_max_compare == None:
    options.z_max_compare = options.z_small - 3
    if options.z_max_compare < 2:
        options.z_max_compare = 2

options.sRNA_features = options.sRNA_features.split(', ')



################################################################################
################################## Curing ######################################
os.chdir(folders.base_folder)
print('Curing wigs')
os.system('python3 {} {} {} {} {}'.format(folders.curing_script, options.wig_folder, genome_size, options.genome_start, options.genome_end))

#
# ### Wait until curing_script is done
if 'temp.txt' in os.listdir(folders.base_folder):
    os.system('rm temp.txt')
else:
    os.wait()


################################################################################
############################## Calling peaks ##################################

os.chdir(folders.base_folder)
print('Calling peaks')
os.system('python3 {} {} {} {} {}'.format(folders.peak_calling_script, options.z_start, options.z_stop, options.z_small, options.processors))
#
# ### Wait until peak_calling script is done
if len([i for i in os.listdir(folders.peak_folder_0) if i != '.DS_Store']) == len([i for i in os.listdir(options.wig_folder) if i != '.DS_Store']) / 2:
    pass
else:
    os.wait()


################################################################################
############################## Make operon file ##################################

## rename files to match operon file

os.chdir(folders.cured_folder)
file_list = sorted(os.listdir())
n = 0
for i in range(0, len(os.listdir()), 2):
    os.system('mv {}{} {}{}_fwd.csv'.format(folders.cured_folder, file_list[i], folders.cured_folder, list_of_conditions[n]))
    os.system('mv {}{} {}{}_rev.csv'.format(folders.cured_folder, file_list[i + 1], folders.cured_folder, list_of_conditions[n]))
    n += 1

os.chdir(folders.peak_folder_0)
file_list = sorted(os.listdir())
n = 0
for i in range(0, len(os.listdir()), 2):
    os.system('mv {}{} {}{}_fwd.csv'.format(folders.peak_folder_0, file_list[i], folders.peak_folder_0, list_of_conditions[n]))
    os.system('mv {}{} {}{}_rev.csv'.format(folders.peak_folder_0, file_list[i + 1], folders.peak_folder_0, list_of_conditions[n]))
    n += 1

list_of_conditions = sorted(list_of_conditions)

os.chdir(folders.base_folder)
if options.operons == None:
    print('Making RendSeq operon file')
    os.system('python3 {} {} {} {}'.format(folders.operon_annotation_script, options.operon_guide, genome_name, options.operon_cutoff))

    # ### Wait until operon script is done
    if len([i for i in os.listdir(folders.operon_folder) if i != '.DS_Store']) == 2:
        pass
    else:
        os.wait()

    os.chdir(folders.operon_folder)

    ### Merge forward and reverse operon files
    both = pd.DataFrame()
    for file in os.listdir():
        if 'GFF3_operon_RendSeq' in file:
            name = file
            df = pd.read_csv(file, sep = '\t', header = 0)
            both = pd.concat([both, df])
            os.system('rm {}'.format(file))
    both.to_csv(name + '_FwdandRev.gff3', sep = '\t', index = False)
else:
    pass


################################################################################
############################## Draw GFF3 ##################################


os.chdir(folders.base_folder)
print('Making GFF3 annotation file')


os.system('python3 {} {} {} {} {} {} {}'.format(folders.assign_and_GFF3_script, folders.peak_folder_0, folders.GFF3_folder_0, options.alternative_operons, genome_name, options.operons, options.processors))



# ### Wait until GFF3 script is done
if len([i for i in os.listdir(folders.GFF3_folder_0) if i != '.DS_Store']) == len([i for i in os.listdir(options.wig_folder) if i != '.DS_Store']) / 2:
    pass
else:
    os.wait()

### Merge forward and reverse annotations and pool conditions

merge_and_pool.merge_and_pool(folders.GFF3_folder_0, folders.GFF3_folder_0, list_of_conditions)


################################################################################
########################## Reevaluate unique peaks ##############################

file = folders.GFF3_folder_0 + 'All_conditions.gff3'
if options.unique != None:
    print('Reevaluating unique peaks')
    start_fwd = options_unique.unique(file, 'start', '+', options.unique)
    start_fwd.unique_and_multiple()

    start_rev = options_unique.unique(file, 'end', '-', options.unique)
    start_rev.unique_and_multiple()

    stop_fwd = options_unique.unique(file, 'end', '+', options.unique)
    stop_fwd.unique_and_multiple()

    stop_rev = options_unique.unique(file, 'start', '-', options.unique)
    stop_rev.unique_and_multiple()

    ## Are there any bidirectional promoters or terminators?
    options_unique.bidirectional(start_fwd.total_list, start_rev.total_list)
    if options_unique.answer == 'yes':
        print('beware of bidirectional promoters')

    options_unique.bidirectional(stop_fwd.total_list, stop_rev.total_list)
    if options_unique.answer == 'yes':
        print('beware of bidirectional terminators')

    unique_starts = options_unique.saver(start_fwd.unique_list, start_rev.unique_list, folders.peak_folder_0, 'unique_starts')
    unique_starts.to_disk()
    unique_stops = options_unique.saver(stop_fwd.unique_list, stop_rev.unique_list, folders.peak_folder_0, 'unique_stops')
    unique_stops.to_disk()

    multiple_starts = options_unique.saver(start_fwd.multiple_list, start_rev.multiple_list, folders.peak_folder_0, 'multiple_starts')
    multiple_starts.to_disk()
    multiple_stops = options_unique.saver(stop_fwd.multiple_list, stop_rev.multiple_list, folders.peak_folder_0, 'multiple_stops')
    multiple_stops.to_disk()

    inst = options_unique.reevaluate(unique_starts.path, unique_stops.path, multiple_starts.path, multiple_stops.path, folders.peak_folder_0, options.z_ii, folders.peak_folder_1)
    inst.func()




################################################################################
############################## Draw GFF3 ##################################

    os.chdir(folders.base_folder)
    # ### Assign to operons and draw GFF3
    print('Making GFF3 annotation file')
    os.system('python3 {} {} {} {} {} {} {}'.format(folders.assign_and_GFF3_script, folders.peak_folder_1, folders.GFF3_folder_1, options.alternative_operons, genome_name, options.operons, options.processors))
    #
    #
    # #
    #
    if len([i for i in os.listdir(folders.GFF3_folder_1) if i != '.DS_Store']) == len([i for i in os.listdir(options.wig_folder) if i != '.DS_Store']) / 2:
        pass
    else:
        os.wait()

    ### Merge forward and reverse annotations and pool conditions
    merge_and_pool.merge_and_pool(folders.GFF3_folder_1, folders.GFF3_folder_1, list_of_conditions)


    #############################################################################################
    ############ Pool annotations where start and stop sites are within 5 bp #####################

    file = folders.GFF3_folder_1 + '/All_conditions.gff3'
    pool_within.pool(file, options.pool, options.wig_folder, folders.GFF3_folder_1, genome_name)
else:
    os.system('mv {} {}'.format(folders.GFF3_folder_0 + '/All_conditions.gff3', folders.GFF3_folder_1 + '/All_conditions.gff3'))
    for file in glob.glob(folders.peak_folder_0 + '*.csv'):
        os.system('mv {} {}/{}'.format(file, folders.peak_folder_1, os.path.basename(file)))


################################################################################
############################ Reevaluate with FDR ################################

file = folders.GFF3_folder_1 + '/All_conditions.gff3'

if options.fdr != 0:
    print('Reevaluating peaks with FDR')
    # fdr.fdr(file, options.fdr, folders.peak_folder_1, folders.peak_folder_2, folders.peak_folder_2, options.z_ii, folders.FDR_tuples_folder)
    os.system('python3 {} {} {} {} {}'.format(folders.fdr_script, folders.peak_folder_1, options.fdr, folders.peak_folder_2, options.processors))
    # for file in glob.glob(folders.peak_folder_1 + '*.csv'):
    #     fdr.fdr_local(file, options.fdr, folders.peak_folder_2)


################################################################################
############################## Draw GFF3 ##################################

    os.chdir(folders.base_folder)
    # ### Assign to operons and draw GFF3
    print('Making GFF3 annotation file')
    os.system('python3 {} {} {} {} {} {} {}'.format(folders.assign_and_GFF3_script, folders.peak_folder_2, folders.GFF3_folder_2, options.alternative_operons, genome_name, options.operons, options.processors))
    #
    #
    # #
    #
    if len([i for i in os.listdir(folders.GFF3_folder_2) if i != '.DS_Store']) == len([i for i in os.listdir(options.wig_folder) if i != '.DS_Store']) / 2:
        pass
    else:
        os.wait()

    ### Merge forward and reverse annotations and pool conditions
    merge_and_pool.merge_and_pool(folders.GFF3_folder_2, folders.GFF3_folder_2, list_of_conditions)


    #############################################################################################
    ############ Pool annotations where start and stop sites are within 5 bp #####################

    # file = folders.GFF3_folder_2 + '/All_conditions.gff3'
    # pool_within.pool(file, options.pool, options.wig_folder, folders.GFF3_folder_2, genome_name)
else:
    os.system('mv {} {}'.format(folders.GFF3_folder_1 + '/All_conditions.gff3', folders.GFF3_folder_2 + '/All_conditions.gff3'))
    for file in glob.glob(folders.peak_folder_1 + '*.csv'):
        os.system('mv {} {}/{}'.format(file, folders.peak_folder_2, os.path.basename(file)))


################################################################################
###################### Reevaluate unique peaks again ##########################

file = folders.GFF3_folder_2 + 'All_conditions.gff3'
if options.unique != None:
    print('Reevaluating unique peaks')
    start_fwd = options_unique.unique(file, 'start', '+', options.unique)
    start_fwd.unique_and_multiple()

    start_rev = options_unique.unique(file, 'end', '-', options.unique)
    start_rev.unique_and_multiple()

    stop_fwd = options_unique.unique(file, 'end', '+', options.unique)
    stop_fwd.unique_and_multiple()

    stop_rev = options_unique.unique(file, 'start', '-', options.unique)
    stop_rev.unique_and_multiple()

    ## Are there any bidirectional promoters or terminators?
    options_unique.bidirectional(start_fwd.total_list, start_rev.total_list)
    if options_unique.answer == 'yes':
        print('beware of bidirectional promoters')

    options_unique.bidirectional(stop_fwd.total_list, stop_rev.total_list)
    if options_unique.answer == 'yes':
        print('beware of bidirectional terminators')

    unique_starts = options_unique.saver(start_fwd.unique_list, start_rev.unique_list, folders.peak_folder_3, 'unique_starts')
    unique_starts.to_disk()
    unique_stops = options_unique.saver(stop_fwd.unique_list, stop_rev.unique_list, folders.peak_folder_3, 'unique_stops')
    unique_stops.to_disk()

    multiple_starts = options_unique.saver(start_fwd.multiple_list, start_rev.multiple_list, folders.peak_folder_3, 'multiple_starts')
    multiple_starts.to_disk()
    multiple_stops = options_unique.saver(stop_fwd.multiple_list, stop_rev.multiple_list, folders.peak_folder_3, 'multiple_stops')
    multiple_stops.to_disk()

    inst = options_unique.reevaluate(unique_starts.path, unique_stops.path, multiple_starts.path, multiple_stops.path, folders.peak_folder_2, options.z_ii, folders.peak_folder_3)
    inst.func()



###############################################################################
############################# Draw GFF3 ##################################

    os.chdir(folders.base_folder)
    # ### Assign to operons and draw GFF3
    print('Making GFF3 annotation file')
    os.system('python3 {} {} {} {} {} {} {}'.format(folders.assign_and_GFF3_script, folders.peak_folder_3, folders.GFF3_folder_3, options.alternative_operons, genome_name, options.operons, options.processors))
    #
    #
    # #
    #
    if len([i for i in os.listdir(folders.GFF3_folder_3) if i != '.DS_Store']) == len([i for i in os.listdir(options.wig_folder) if i != '.DS_Store']) / 2:
        pass
    else:
        os.wait()

    ### Merge forward and reverse annotations and pool conditions
    merge_and_pool.merge_and_pool(folders.GFF3_folder_3, folders.GFF3_folder_3, list_of_conditions)


    #############################################################################################
    ############ Pool annotations where start and stop sites are within 5 bp #####################

    file = folders.GFF3_folder_3 + '/All_conditions.gff3'
    pool_within.pool(file, options.pool, options.wig_folder, folders.GFF3_folder_3, genome_name)
else:
    os.system('mv {} {}'.format(folders.GFF3_folder_2 + '/All_conditions.gff3', folders.GFF3_folder_3 + '/All_conditions.gff3'))
    for file in glob.glob(folders.peak_folder_2 + '*.csv'):
        os.system('mv {} {}/{}'.format(file, folders.peak_folder_3, os.path.basename(file)))


    ###########################################################################################
    ########## Supply using alternative RNA-seq methods (e.g. cap-seq and term-seq) ###########

### Depending on options choosen move GFF3 files to GFF3_folder_3 for further processing
if options.fdr == 0 and options.unique == None:
    for file in glob.glob(folders.GFF3_folder_1 + '*.gff3'):
        os.system('mv {} {}{}'.format(file, folders.GFF3_folder_3, os.path.basename(file)))

if options.fdr != 0 and options.unique == None:
    for file in glob.glob(folders.GFF3_folder_2 + '*.gff3'):
        os.system('mv {} {}{}'.format(file, folders.GFF3_folder_3, os.path.basename(file)))


cap_term = {}
if options.cap_seq != None:
    cap = pd.read_csv(options.cap_seq, sep = ',', header = 0)
    cap_fwd = cap['location'][cap['strand'] == '+'].tolist()
    cap_rev = cap['location'][cap['strand'] == '-'].tolist()
    cap_term['starts_fwd'] = cap_fwd
    cap_term['starts_rev'] = cap_rev
if options.term_seq != None:
    term = pd.read_csv(options.term_seq, sep = ',', header = 0)
    if 'strand' in term.columns:
        term_fwd = term['location'][term['strand'] == '+'].tolist()
        term_rev = term['location'][term['strand'] == '-'].tolist()
        cap_term['stops_fwd'] = term_fwd
        cap_term['stops_rev'] = term_rev
    else:
        term_all = term['location'].tolist()
        cap_term['stops_all'] = term_all

operon_starts_stops = {}
if options.operon_starts != None:
    cap = pd.read_csv(options.operon_starts, sep = ',', header = 0)
    cap_fwd = cap['location'][cap['strand'] == '+'].tolist()
    cap_rev = cap['location'][cap['strand'] == '-'].tolist()
    operon_starts_stops['starts_fwd'] = cap_fwd
    operon_starts_stops['starts_rev'] = cap_rev
if options.operon_stops != None:
    term = pd.read_csv(options.operon_stops, sep = ',', header = 0)
    term_fwd = term['location'][term['strand'] == '+'].tolist()
    term_rev = term['location'][term['strand'] == '-'].tolist()
    operon_starts_stops['stops_fwd'] = term_fwd
    operon_starts_stops['stops_rev'] = term_rev

if options.cap_seq != None or options.term_seq != None:
    if options.old_annotation != None:
        print('Searching for peaks using alternative peak files')
        np.save(os.path.join(folders.analysis_folder, 'cap_term.npy'), cap_term)
        for file in glob.glob(folders.peak_folder_3 + '*.csv'):
            condition = os.path.basename(file)[:-8]
            gff_file = folders.GFF3_folder_3 + condition + '.gff3'
            lone_peaks = Draw_GFF3_supply_w_capterm.missing_peaks(file, gff_file, folders.analysis_folder + 'cap_term.npy', genome_name, condition, options.old_annotation, options.search_distance, folders.GFF3_folder_3, 1)
            lone_peaks.supply_lone_peaks()
            lone_peaks = Draw_GFF3_supply_w_capterm.missing_peaks(file, gff_file, folders.analysis_folder + 'cap_term.npy', genome_name, condition, options.old_annotation, options.search_distance, folders.GFF3_folder_3, 2)
            lone_peaks.supply_remaining_peaks()

if options.operon_starts != None or options.operon_stops != None:
    if options.old_annotation != None:
        print('Searching for peaks using operon start and stop sites')
        np.save(os.path.join(folders.analysis_folder, 'operon_starts_stops.npy'), operon_starts_stops)
        for file in glob.glob(folders.peak_folder_3 + '*.csv'):
            condition = os.path.basename(file)[:-8]
            gff_file = folders.GFF3_folder_3 + condition + '.gff3'
            lone_peaks = Draw_GFF3_supply_w_capterm.missing_peaks(file, gff_file, folders.analysis_folder + 'operon_starts_stops.npy', genome_name, condition, options.old_annotation, options.search_distance_0, folders.GFF3_folder_3, 2)
            lone_peaks.supply_lone_peaks()
            lone_peaks = Draw_GFF3_supply_w_capterm.missing_peaks(file, gff_file, folders.analysis_folder + 'operon_starts_stops.npy', genome_name, condition, options.old_annotation, options.search_distance_0, folders.GFF3_folder_3, 2)
            lone_peaks.supply_remaining_peaks()

if options.cap_seq != None or options.term_seq != None:
    if options.old_annotation != None:
        for file in glob.glob(folders.peak_folder_3 + '*.csv'):
            condition = os.path.basename(file)[:-8]
            gff_file = folders.GFF3_folder_3 + condition + '.gff3'
            lone_peaks = Draw_GFF3_supply_w_capterm.missing_peaks(file, gff_file, folders.analysis_folder + 'cap_term.npy', genome_name, condition, options.old_annotation, options.search_distance, folders.GFF3_folder_3, 2)
            lone_peaks.supply_half_peaks()

if options.operon_starts != None or options.operon_stops != None:
    if options.old_annotation != None:
        for file in glob.glob(folders.peak_folder_3 + '*.csv'):
            condition = os.path.basename(file)[:-8]
            gff_file = folders.GFF3_folder_3 + condition + '.gff3'
            lone_peaks = Draw_GFF3_supply_w_capterm.missing_peaks(file, gff_file, folders.analysis_folder + 'operon_starts_stops.npy', genome_name, condition, options.old_annotation, options.search_distance_0, folders.GFF3_folder_3, 2)
            lone_peaks.supply_half_peaks()

if options.operon_starts != None or options.operon_stops != None:
    if options.old_annotation != None:
        for file in glob.glob(folders.peak_folder_0 + '*.csv'):
            condition = os.path.basename(file)[:-8]
            gff_file = folders.GFF3_folder_3 + condition + '.gff3'
            lone_peaks = Draw_GFF3_supply_w_capterm.missing_peaks(file, gff_file, folders.analysis_folder + 'operon_starts_stops.npy', genome_name, condition, options.old_annotation, options.search_distance_0, folders.GFF3_folder_3, 2)
            lone_peaks.supply_half_peaks_operons()


os.system('python3 {} {} {} {}'.format(folders.delete_redundant_script, folders.GFF3_folder_3, folders.GFF3_folder_4, options.processors))

# ### Wait until GFF3 script is done
if len([i for i in os.listdir(folders.GFF3_folder_4) if i != '.DS_Store']) -1 == len([i for i in os.listdir(options.wig_folder) if i != '.DS_Store']) / 4:
    pass
else:
    os.wait()

merge_and_pool.pool(folders.GFF3_folder_4)
file = folders.GFF3_folder_4 + '/All_conditions.gff3'
pool_within.pool(file, options.pool, options.wig_folder, folders.GFF3_folder_4, genome_name)

#

# ##############################################################################
# ##############################################################################
# ################# Transfer features from Nicolas et al.#######################
# ## The different features in the Nicolas et al. annotation are
# ## ['CDS', 'transcript', 'sequence_feature', 'tRNA', 'rRNA', 'tmRNA', 'sRNA']
# ## We transfer annotations to the Rend-seq file if there is an overlap with a
# ## Nicolas et al. annotation.
# ## If a Nicolas CDS is found completely within a Rend-seq transcript then that Rend-seq transcript
# ## will be assigned the feature CDS.
# ## For all Rend-seq transcripts that are not CDS, they will be given the same feature
# ## as any nicolas gene with which they overlap.
# ## If a Rend-seq trancript doesn't overlap with any Nicolas genes it will be designated a novel transcript
# ## If a non-CDS gene overlaps with multiple Nicolas genes if will be given the feature 'transcript'

### Distance between Rend-seq and Nicolas et al. ###################
# # For all the genes that overlap with a Rend-seq transcript determine the
# # minimum distance to start and and stop
#



if options.old_annotation != None:
    for file in glob.glob(folders.GFF3_folder_4 + '*.gff3'):
        post_processing.scan_halfs(file)
    post = post_processing.compare(options.old_annotation, options.source_exclude, folders.GFF3_folder_4, options.genome_start, options.genome_end, genome_name, folders.analysis_folder, old_annotation_name, options.sRNA_features, options.stat_cond)
    post.phase()
    print('Transfering features and names from old annotation file')
    post.features()
    print('Categorizing sRNAs')
    post.categorize_sRNAs()
    print('Categorizing operons')
    post.categorize_operons()
    print('Finding genes in old annotation file not covered by Rend-seq')
    post.not_covered()
    print('Making new annotation supplied with old annotation where there is no coverage')
    post.supply_cat_operon()
    post.supply_raw()
    for file in glob.glob(folders.GFF3_folder_4 + '*.gff3'):
        post_processing.delete_1_bp(file)
        post_processing.relabel_novel(file, options.old_annotation, 0)
    post_processing.relabel_novel(f'{folders.analysis_folder}/Statistics_{options.stat_cond}.txt', options.old_annotation, 1)


# Figure 5 shows a screenshot of JBrowse that shows an example of what the annotation file produced by this
# pipeline contains.
#
#
# ### move everything to provided output_folder (--output_folder)
try:
    os.chdir(options.output_folder)
    if '.DS_Store' in os.listdir():
        os.remove('.DS_Store')
    if 'Terminal_input.py' in os.listdir():
        os.remove('Terminal_input.py')
    for folder in os.listdir():
        shutil.rmtree(folder, ignore_errors=False, onerror=None)
except:
    pass

### remove parent= in annotations as it don't function in Jbrowse
if options.old_annotation != None:
    for file in glob.glob(folders.GFF3_folder_4 + '/*.gff3'):
        df = pd.read_csv(file, sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
        attributes = df['attributes'].tolist()
        new_attributes = []
        for item in attributes:
            alist = item.split(';')
            newlist = []
            for item in alist:
                if item.find('Parent=') == 0:
                    pass
                else:
                    newlist.append(item)
            astring = ';'.join(newlist)
            new_attributes.append(astring)
        df['attributes'] = new_attributes
        df.to_csv(os.path.join(folders.GFF3_folder_4, os.path.basename(file)), sep = '\t', header = None, index = False)


try:
    os.makedirs(options.output_folder)
except:
    pass



if options.statistics == True:
    print('Making statistics figures')
    statistical_analysis.read_density(folders.cured_folder, list_of_conditions, folders.analysis_folder)
    statistical_analysis.plot_n_peaks([i for i in glob.glob(folders.GFF3_folder_4 + '*.gff3') if 'all_rend_peaks' not in i and 'categorized_operons' not in i and 'supplied_with_' not in i], folders.analysis_folder)
    if options.operon_starts != None and options.operon_stops != None:
        statistical_analysis.plot_coverage(options.operon_starts, options.operon_stops, folders.GFF3_folder_4 + 'All_conditions.gff3', options.search_distance_0, folders.analysis_folder)
    if options.cap_seq != None and options.term_seq != None:
        try:
            statistical_analysis.plot_agreement_w_CapTerm(options.cap_seq, options.term_seq, folders.GFF3_folder_4 + '/' + options.compare_condition + '.gff3', folders.analysis_folder)
        except:
            print('Could not run part of statistical analysis to plot agreement with alternative coordinates')
    if options.feature_distance != None:
        try:
            statistical_analysis.distance(folders.analysis_folder + '/Statistics_All_conditions.txt', folders.distance_folder, options.feature_distance)
        except:
            print('Could not run part of statistical analysis to plot distance to alternative coordinates. Make sure your --feature_distance file is correct.')
    if options.compare_peaks_across_samples != None:
        try:
            statistical_analysis.differentially_expressed_peaks(folders.GFF3_folder_4, options.compare_peaks_across_samples, folders.peak_folder_0, options.z_min_compare, options.z_max_compare, options.pool, options.density_compare, folders.differentially_expressed_genes_Overlap) #, folders.differentially_expressed_genes_All, folders.differentially_expressed_genes_HighConfidence, folders.differentially_expressed_genes_HighDensity, folders.differentially_expressed_genes_Overlap_lessStringent)
        except:
            print('Could not find differentially expressed genes. Make sure your --compare_peaks_across_samples file is correct.')



if options.output == 'All':
    os.system('mv {} {}'.format(folders.pipeline_folder, options.output_folder))
else:
    os.system('mv {} {}'.format(folders.GFF3_folder_4, options.output_folder))
    try:
        shutil.rmtree(folders.pipeline_folder)
    except:
        pass

os.chdir(folders.base_folder)
if options.output == 'GTF':
    for file in glob.glob(options.output_folder + '/Annotation_files/*.gff3'):
        os.system('python3 {} -f {} -o {}{}'.format(folders.GFF_to_GTF_script, file, file[:-4], 'gtf'))


try:
    shutil.copy(folders.base_folder + '/EXAMPLE_Terminal_input.py', options.output_folder + '/EXAMPLE_Terminal_input.py')
except:
    pass
