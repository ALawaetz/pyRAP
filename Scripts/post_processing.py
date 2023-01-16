import os
import pandas as pd
import glob
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
#	post_processing.py
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


# ##############################################################################
# ##############################################################################
# ################# Transfer features from old_annotation file #######################
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
class compare:
    def __init__(self, old_annotation, rend_folder, genome_start, genome_end, genome_name, analysis_folder, old_annotation_name):
        self.df_Nicolas = pd.read_csv(old_annotation, sep = '\t', header = 0)
        self.rend_folder = rend_folder
        self.genome_start = genome_start
        self.genome_end = genome_end
        self.genome_name = genome_name
        self.analysis_folder = analysis_folder
        self.old_annotation_name = old_annotation_name


    def phase(self):
        for file in glob.glob(self.rend_folder + '*.gff3'):
            df_rend = pd.read_csv(file, sep = '\t', header = 0)
            new_df = pd.DataFrame()
            for index, row in df_rend.iterrows():
                if row['phase'] == '.':
                    new_df = pd.concat([new_df, df_rend[df_rend.index == index]])
                else:
                    mydata = {'seqID': [row['seqID']], 'source': [row['source']], 'feature': [row['feature']], 'start': [row['start']], 'end': [row['end']], 'score': [row['score']], 'strand': [row['strand']], 'phase': ['dot'], 'attributes': ['{};Note={}'.format(row['attributes'], row['phase'])]}
                    df_1 = pd.DataFrame(data = mydata)
                    new_df = pd.concat([new_df, df_1])
            new_df.to_csv(os.path.join(self.rend_folder, os.path.basename(file)), sep = '\t', index = False)

    def features(self):

        self.df_Nicolas = self.df_Nicolas[(self.df_Nicolas['start'] >= self.genome_start) & (self.df_Nicolas['end'] <= self.genome_end)]
        self.df_Nicolas = self.df_Nicolas.sort_values(['start', 'end'], axis = 0, ascending = True)
        self.df_Nicolas.reset_index(inplace = True, drop = True)

        for file in glob.glob(self.rend_folder + '*.gff3'):
            Rendseq_file = file
            try:
                df_rend = pd.read_csv(Rendseq_file, sep = '\t', header = 0)
                df_rend = df_rend.sort_values(['start', 'end'], axis = 0, ascending = True)
                df_rend.reset_index(inplace = True, drop = True)
            except:
                print('''No annotations made for {}
                Try lowering your cutoff values using the z-flags'''.format(os.path.basename(file)))
                df_rend = pd.DataFrame()


            rend_feature = pd.DataFrame()

            start_col = []
            stop_col = []
            for index, row in df_rend.iterrows():
                R_strand = row['strand']
                R_start = row['start']
                R_end = row['end']
                ### If there is an overlap then one of these conditions will provide a dataframe in df_rend with a length above 0
                ### This method can be discussed as a RendSeq annotation might overlap with a Nicolas et al. annotation even though they are not the same genes
                ### e.g. yabG gene
                conditions = (R_start <= self.df_Nicolas['start']) & (R_end > self.df_Nicolas['start']) | (R_start >= self.df_Nicolas['start']) & (R_start < self.df_Nicolas['end']) | (R_start == self.df_Nicolas['start']) & (R_end == self.df_Nicolas['end'])
                if len(self.df_Nicolas[conditions][self.df_Nicolas[conditions]['strand'] == R_strand]) > 0:
                    ### Determine shortest distance to start or stop sites in Nicolas et al.
                    cond_start = []
                    cond_start_abs = []
                    cond_stop = []
                    cond_stop_abs = []
                    start_pos = self.df_Nicolas[conditions]['start'][self.df_Nicolas[conditions]['strand'] == R_strand].tolist()
                    end_pos = self.df_Nicolas[conditions]['end'][self.df_Nicolas[conditions]['strand'] == R_strand].tolist()
                    for item in start_pos:
                        cond_start.append(R_start - item)
                        cond_start_abs.append(abs(R_start - item))
                    for item in end_pos:
                        cond_stop.append(R_end - item)
                        cond_stop_abs.append(abs(R_end - item))
                    if  R_strand == '+':
                        min_start = cond_start[cond_start_abs.index(min(cond_start_abs))]
                        min_stop = cond_stop[cond_stop_abs.index(min(cond_stop_abs))]
                    else:
                        min_start = -cond_stop[cond_stop_abs.index(min(cond_stop_abs))]
                        min_stop = -cond_start[cond_start_abs.index(min(cond_start_abs))]
                    start_col.append(min_start)
                    stop_col.append(min_stop)

                    ### make a list of the names of the genes that overlap
                    genes = ''
                    features = []
                    gene_count = 0
                    for item in self.df_Nicolas[conditions]['attributes'][self.df_Nicolas[conditions]['strand'] == R_strand].tolist():
                        alist = item.split(';')
                        z = 0
                        for a in alist:
                            if a.find('Name=') == 0:
                                if gene_count == 0:
                                    genes += a[5:]
                                    gene_count += 1
                                    z += 1
                                else:
                                    genes += '_{}'.format(a[5:])
                                    z += 1
                            else:
                                pass
                        if z == 0:
                            for a in alist:
                                if a.find('ID=') == 0:
                                    if gene_count == 0:
                                        genes += a[3:]
                                        gene_count += 1
                                        z += 1
                                    else:
                                        genes += '_{}'.format(a[3:])
                                        z += 1
                                else:
                                    pass
                        if z == 0:
                            for a in alist:
                                if a.find('symbol=') == 0:
                                    if gene_count == 0:
                                        genes += a[7:]
                                        gene_count += 1
                                        z += 1
                                    else:
                                        genes += '_{}'.format(a[7:])
                                        z += 1
                                else:
                                    pass
                        if z == 0:
                            for a in alist:
                                if a.find('Description=') == 0:
                                    if gene_count == 0:
                                        genes += a[12:]
                                        gene_count += 1
                                        z += 1
                                    else:
                                        genes += '_{}'.format(a[12:])
                                        z += 1
                                else:
                                    pass
                        if z == 0:
                            for a in alist:
                                if a.find('product=') == 0:
                                    if gene_count == 0:
                                        genes += a[8:]
                                        gene_count += 1
                                        z += 1
                                    else:
                                        genes += '_{}'.format(a[8:])
                                        z += 1
                                else:
                                    pass

                    if len(genes) == 0:
                        genes = ''
                        for item in self.df_Nicolas[conditions]['feature'][self.df_Nicolas[conditions]['strand'] == R_strand].tolist():
                            if gene_count == 0:
                                genes += item
                                gene_count += 1
                            else:
                                genes += '_{}'.format(item)



                    for item in self.df_Nicolas[conditions]['feature'][self.df_Nicolas[conditions]['strand'] == R_strand].tolist():
                        features.append(item)
                    ### For coding sequences there has to be a complete overlap or more
                    conditions_CDS = (R_start <= self.df_Nicolas['start']) & (R_end >= self.df_Nicolas['end']) & (R_strand == self.df_Nicolas['strand']) | (R_start == self.df_Nicolas['start']) & (R_end == self.df_Nicolas['end']) & (R_strand == self.df_Nicolas['strand'])
                    if row['phase'] == 'dot':
                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ['name={}; {}; start_dist={}; stop_dist={}; Type=Only one peak resolved'.format(genes, row['attributes'], min_start, min_stop)]}
                        df_1 = pd.DataFrame(data = mydata)
                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                    elif 'CDS' in self.df_Nicolas[conditions_CDS]['feature'].tolist():
                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['CDS'], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ['name={}; {}; start_dist={}; stop_dist={}'.format(genes, row['attributes'], min_start, min_stop)]}
                        df_1 = pd.DataFrame(data = mydata)
                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                    else:
                        ### There might still be a CDS with a partial overlap
                        if len(list(set(features))) == 2 and 'CDS' in features:
                            thefeature = [i for i in features if i != 'CDS'][0]
                            mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [thefeature], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ['name={}; {}; start_dist={}; stop_dist={}'.format(genes, row['attributes'], min_start, min_stop)]}
                            df_1 = pd.DataFrame(data = mydata)
                            rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                        elif len(list(set(features))) > 2:
                            ### more than two gene feature overlap and we will just write transcipt
                            ### alternatively we could pick the feature of the gene with greatest overlap
                            mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['transcript'], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ['name={}; {}; start_dist={}; stop_dist={}'.format(genes, row['attributes'], min_start, min_stop)]}
                            df_1 = pd.DataFrame(data = mydata)
                            rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                        else:
                            if features[0] != 'CDS':
                                mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['{}'.format(features[0])], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ['name={}; {}; start_dist={}; stop_dist={}'.format(genes, row['attributes'], min_start, min_stop)]}
                                df_1 = pd.DataFrame(data = mydata)
                                rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                            else:
                                ### else the overlap is with a CDS but the Rend-seq transcript doesn't cover the entire CDS so we designate it putative_sRNA
                                ### unless the transcript covers the 3' end of a Nicolas et al annotation in which case it can still potentially encode an ORF in the
                                ### same reading frame so we will designate it CDS and write in attributes 'alternative start site'
                                if R_strand == '+':
                                    for end in self.df_Nicolas[conditions]['end'][self.df_Nicolas[conditions]['strand'] == R_strand].tolist():
                                        count = 0
                                        if end <= R_end:
                                            mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['CDS'], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ['name={}; {}; start_dist={}; stop_dist={}; Type=Alternative start site'.format(genes, row['attributes'], min_start, min_stop)]}
                                            df_1 = pd.DataFrame(data = mydata)
                                            rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                            count += 1
                                            break
                                        else:
                                            pass
                                    if count == 0:
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['putative_sRNA'], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ['name={}; {}; start_dist={}; stop_dist={}'.format(genes, row['attributes'], min_start, min_stop)]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    else:
                                        pass
                                else:
                                    for start in self.df_Nicolas[conditions]['start'][self.df_Nicolas[conditions]['strand'] == R_strand].tolist():
                                        count = 0
                                        if start >= R_start:
                                            mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['CDS'], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ['name={}; {}; start_dist={}; stop_dist={}; Type=Alternative start site'.format(genes, row['attributes'], min_start, min_stop)]}
                                            df_1 = pd.DataFrame(data = mydata)
                                            rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                            count += 1
                                            break
                                        else:
                                            pass
                                    if count == 0:
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['putative_sRNA'], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ['name={}; {}; start_dist={}; stop_dist={}'.format(genes, row['attributes'], min_start, min_stop)]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    else:
                                        pass
                else:
                    ### Novel transcript
                    ### The Nicolas at al. annoations has no transcipts assigned the feature 'gene'. Therefore we assign all new transcripts the feature 'gene' whereby
                    ### we can easily count the number of new transcripts.
                    start_col.append(0)
                    stop_col.append(0)
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': ['Novel_transcript'], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ['name=Novel_transcript; {}'.format(row['attributes'])]}
                    df_1 = pd.DataFrame(data = mydata)
                    rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)

            rend_feature.to_csv(os.path.join(self.rend_folder, os.path.basename(file)), sep = '\t', index = False)


            ## Make a statistics file
            if os.path.basename(file) == 'All_conditions.gff3':
                rend_feature['start_dist'] = start_col
                rend_feature['stop_dist'] = stop_col
                rend_feature.to_csv(os.path.join(self.analysis_folder, 'Statistics_All_conditions.txt'), sep = '\t', index = False)
    #
# ############################################################################
# ############################################################################
# ########################## Categorize sRNAs ################################
# # Sub-categorize sRNAs and putative sRNAs into derived from 5'UTRs, 3'UTRs, intragenic, intergenic,
# # and independent
### Rend-seq Transcripts that are not CDS but overlap with a Rend-seq transcript will be categorized as either
# # 5'UTRs, 3'UTRs, intragenic, intergenic or intragenic/intergenic sRNAs
# # For the Nicolas genes which overlap with the Rend-seq CDS that overlap with the sRNA in question, if start of the sRNA is upstream
# # the Nicolas CDS start, then it is 5'
# # If the sRNA stop is downstream the Nicolas CDS stop, then it is 3'
# # if it is within the CDS it is intragenic
# # If it is between two Nicolas CDS, then it is intergenic
# # If the sRNA doesn't overlap with a Rend-seq CDS, then it is independent.

    def categorize_sRNAs(self):
        self.df_Nicolas = self.df_Nicolas.sort_values(['start', 'end'], axis = 0, ascending = True)
        self.df_Nicolas.reset_index(inplace = True, drop = True)
        for file in glob.glob(self.rend_folder + '*.gff3'):
            try:
                df_rend = pd.read_csv(file, sep = '\t', header = 0)
                df_rend = df_rend.sort_values(['start', 'end'], axis = 0, ascending = True)
                df_rend.reset_index(inplace = True, drop = True)
            except:
                print('''No annotaions made for {}
                Consider lowering your z-score using one of the z flags or analyse a larger genomic region if using the --genome_start and end flags.
                For help type: python3 1_Rend_seq_pipeline_22August2022.py --help'''.format(os.path.basename(file)))
                continue

            rend_feature = pd.DataFrame()
            sRNA_type = []

            for index, row in df_rend.iterrows():
                if row['feature'] == 'sRNA' or row['feature'] == 'putative_sRNA' or row['feature'] == 'antisense_RNA' or row['feature'] == 'ncRNA' or row['feature'] == 'misc_RNA':
                    R_strand = row['strand']
                    R_start = row['start']
                    R_end = row['end']
                    ### If there is an overlap then one of these conditions will provide a dataframe in df_rend with a length above 0
                    ### Does the Rend-Seq annotation overlap with any other Rend-seq annotations
                    conditions = (R_start <= df_rend['start']) & (R_end > df_rend['start']) | (R_start >= df_rend['start']) & (R_start < df_rend['end']) | (R_start == df_rend['start']) & (R_end == df_rend['end'])
                    conditions_N = (R_start <= self.df_Nicolas['start']) & (R_end > self.df_Nicolas['start']) | (R_start >= self.df_Nicolas['start']) & (R_start < self.df_Nicolas['end']) | (R_start == self.df_Nicolas['start']) & (R_end == self.df_Nicolas['end'])
                    if len(df_rend[conditions][df_rend[conditions]['strand'] == R_strand]) > 1:
                        if 'CDS' in df_rend[conditions]['feature'][df_rend[conditions]['strand'] == R_strand].tolist():
                            ### Being here means an sRNA that is either derived from 5'UTR, 3'UTR, intragenic or intergenic
                            ### because it overlaps with a CDS Rend-seq transcript.
                            ### Find start and end of the overlapping CDS
                            CDS = df_rend[conditions][(df_rend[conditions]['strand'] == R_strand) & (df_rend[conditions]['feature'] == 'CDS')]
                            R_start_over = min(CDS['start'].tolist())
                            R_end_over = max(CDS['end'].tolist())
                            R_strand_over = R_strand
                            conditions = (R_start_over <= self.df_Nicolas['start']) & (R_end_over > self.df_Nicolas['start']) | (R_start_over >= self.df_Nicolas['start']) & (R_start_over < self.df_Nicolas['end']) | (R_start_over == self.df_Nicolas['start']) & (R_end_over == self.df_Nicolas['end'])
                            overlap = self.df_Nicolas[conditions][(self.df_Nicolas[conditions]['strand'] == R_strand_over) & (self.df_Nicolas[conditions]['feature'] == 'CDS')]
                            ### now see if the row of the iteration overlaps with any CDS in the overlap
                            conditions = (R_start <= overlap['start']) & (R_end > overlap['start']) | (R_start >= overlap['start']) & (R_start < overlap['end']) | (R_start == overlap['start']) & (R_end == overlap['end'])
                            if len(overlap[conditions][overlap[conditions]['strand'] == R_strand]) == 0:
                                    ## The sRNA doesn't overlap with any Nicolas CDS
                                ### sRNA is either 5', or 3' or intergenic
                                if R_strand == '+':
                                    if R_start < min(overlap['start']):
                                        ### 5'UTR
                                        sRNA_type.append("5'UTR")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=5'UTR".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    elif R_end > max(overlap['end']):
                                        ### 3'UTR
                                        sRNA_type.append("3'UTR")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=3'UTR".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    else:
                                        ### intergenic
                                        sRNA_type.append("intergenic")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=intergenic".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                else:
                                    if R_end > max(overlap['end']):
                                        ### 5'UTR
                                        sRNA_type.append("5'UTR")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=5'UTR".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    elif R_start < min(overlap['start']):
                                        ### 3'UTR
                                        sRNA_type.append("3'UTR")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=3'UTR".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    else:
                                        ### intergenic
                                        sRNA_type.append("intergenic")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=intergenic".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                            elif len(overlap[conditions][overlap[conditions]['strand'] == R_strand]) == 1:
                                ### 5', 3' or intragenic
                                if R_strand == '+':
                                    if R_start < min(overlap[conditions]['start'][overlap[conditions]['strand'] == R_strand]):
                                        ### 5'UTR
                                        sRNA_type.append("5'UTR")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=5'UTR".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    elif R_end > max(overlap[conditions]['end'][overlap[conditions]['strand'] == R_strand]):
                                        ### 3'UTR
                                        sRNA_type.append("3'UTR")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=3'UTR".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    else:
                                        ### intragenic
                                        sRNA_type.append("intragenic")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=intragenic".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                else:
                                    if R_end > max(overlap[conditions]['end'][overlap[conditions]['strand'] == R_strand]):
                                        ### 5'UTR
                                        sRNA_type.append("5'UTR")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=5'UTR".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    elif R_start < min(overlap[conditions]['start'][overlap[conditions]['strand'] == R_strand]):
                                        ### 3'UTR
                                        sRNA_type.append("3'UTR")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=3'UTR".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    else:
                                        ### intragenic
                                        sRNA_type.append("intragenic")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=intragenic".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                            else:
                                ### overlap with more than one CDS means inter/intragenic
                                sRNA_type.append("inter/intragenic")
                                mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=inter/intragenic".format(row['attributes'])]}
                                df_1 = pd.DataFrame(data = mydata)
                                rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                        ### if the transcript of the iteration overlaps with other Rend-seq transcript, but none of these are CDS
                        ### then it is an sRNA with isoforms if more than one of the overlapping transcripts are sRNAs or putative sRNAs
                        ### Can also be sRNA isoforms of protein (see example with rapA)
                        ### Therefore check if annotation overlaps any Nicolas CDS
                        elif len(self.df_Nicolas[conditions_N][(self.df_Nicolas[conditions_N]['strand'] == R_strand) & (self.df_Nicolas[conditions_N]['feature'] == 'CDS')]) == 0:
                                if df_rend[conditions]['feature'][df_rend[conditions]['strand'] == R_strand].tolist().count('sRNA') + df_rend[conditions]['feature'][df_rend[conditions]['strand'] == R_strand].tolist().count('putative_sRNA') + df_rend[conditions]['feature'][df_rend[conditions]['strand'] == R_strand].tolist().count('antisense_RNA') + df_rend[conditions]['feature'][df_rend[conditions]['strand'] == R_strand].tolist().count('ncRNA') + df_rend[conditions]['feature'][df_rend[conditions]['strand'] == R_strand].tolist().count('misc_RNA') + df_rend[conditions]['feature'][df_rend[conditions]['strand'] == R_strand].tolist().count('Novel_transcript')  > 1:
                                    sRNA_type.append("sRNA_isoform")
                                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=independent_w_isoform".format(row['attributes'])]}
                                    df_1 = pd.DataFrame(data = mydata)
                                    rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                else:
                                    ### sRNA that overlaps with genes other than sRNAs and CDS
                                    sRNA_type.append("independent_overlap")
                                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=non_sRNA/non_CDS_overlap".format(row['attributes'])]}
                                    df_1 = pd.DataFrame(data = mydata)
                                    rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                        else:
                            ### the sRNA doesn't overlap with any other RendSeq CDS annotations but might still overlap with a Nicolas annotation
                            ### if it doesn't overlap with any Nicolas CDS annotations then it is independent
                            conditions = (R_start <= self.df_Nicolas['start']) & (R_end > self.df_Nicolas['start']) | (R_start >= self.df_Nicolas['start']) & (R_start < self.df_Nicolas['end']) | (R_start == self.df_Nicolas['start']) & (R_end == self.df_Nicolas['end'])
                            overlap = self.df_Nicolas[conditions][(self.df_Nicolas[conditions]['strand'] == R_strand) & (self.df_Nicolas[conditions]['feature'] == 'CDS')]
                            if len(overlap) == 1:
                                if R_strand == '+':
                                    if R_start < min(overlap['start']):
                                        ### 5'UTR
                                        sRNA_type.append("5'UTR")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=5'UTR".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    elif R_end > max(overlap['end']):
                                        ### 3'UTR
                                        sRNA_type.append("3'UTR")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=3'UTR".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    else:
                                        ### intragenic
                                        sRNA_type.append("intragenic")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=intragenic".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                else:
                                    if R_end > max(overlap['end']):
                                        ### 5'UTR
                                        sRNA_type.append("5'UTR")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=5'UTR".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    elif R_start < min(overlap['start']):
                                        ### 3'UTR
                                        sRNA_type.append("3'UTR")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=3'UTR".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                    else:
                                        ### intragenic
                                        sRNA_type.append("intragenic")
                                        mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=intragenic".format(row['attributes'])]}
                                        df_1 = pd.DataFrame(data = mydata)
                                        rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                            elif len(overlap) > 1:
                                ### overlap with more than one CDS means inter/intragenic
                                sRNA_type.append("inter/intragenic")
                                mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=inter/intragenic".format(row['attributes'])]}
                                df_1 = pd.DataFrame(data = mydata)
                                rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                            else:
                                sRNA_type.append("independent")
                                mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=independent".format(row['attributes'])]}
                                df_1 = pd.DataFrame(data = mydata)
                                rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                    else:
                        ### the sRNA doesn't overlap with any other RendSeq CDS annotations but might still overlap with a Nicolas annotation
                        ### if it doesn't overlap with any Nicolas CDS annotations then it is independent
                        conditions = (R_start <= self.df_Nicolas['start']) & (R_end > self.df_Nicolas['start']) | (R_start >= self.df_Nicolas['start']) & (R_start < self.df_Nicolas['end']) | (R_start == self.df_Nicolas['start']) & (R_end == self.df_Nicolas['end'])
                        overlap = self.df_Nicolas[conditions][(self.df_Nicolas[conditions]['strand'] == R_strand) & (self.df_Nicolas[conditions]['feature'] == 'CDS')]
                        if len(overlap) == 1:
                            if R_strand == '+':
                                if R_start < min(overlap['start']):
                                    ### 5'UTR
                                    sRNA_type.append("5'UTR")
                                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=5'UTR".format(row['attributes'])]}
                                    df_1 = pd.DataFrame(data = mydata)
                                    rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                elif R_end > max(overlap['end']):
                                    ### 3'UTR
                                    sRNA_type.append("3'UTR")
                                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=3'UTR".format(row['attributes'])]}
                                    df_1 = pd.DataFrame(data = mydata)
                                    rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                else:
                                    ### intragenic
                                    sRNA_type.append("intragenic")
                                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=intragenic".format(row['attributes'])]}
                                    df_1 = pd.DataFrame(data = mydata)
                                    rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                            else:
                                if R_end > max(overlap['end']):
                                    ### 5'UTR
                                    sRNA_type.append("5'UTR")
                                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=5'UTR".format(row['attributes'])]}
                                    df_1 = pd.DataFrame(data = mydata)
                                    rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                elif R_start < min(overlap['start']):
                                    ### 3'UTR
                                    sRNA_type.append("3'UTR")
                                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=3'UTR".format(row['attributes'])]}
                                    df_1 = pd.DataFrame(data = mydata)
                                    rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                                else:
                                    ### intragenic
                                    sRNA_type.append("intragenic")
                                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=intragenic".format(row['attributes'])]}
                                    df_1 = pd.DataFrame(data = mydata)
                                    rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                        elif len(overlap) > 1:
                            ### overlap with more than one CDS means inter/intragenic
                            sRNA_type.append("inter/intragenic")
                            mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=inter/intragenic".format(row['attributes'])]}
                            df_1 = pd.DataFrame(data = mydata)
                            rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                        else:
                            sRNA_type.append("independent")
                            mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [R_start], 'end': [R_end], 'score': [row['score']], 'strand': [R_strand], 'phase': [row['phase']], 'attributes': ["{}; sRNA_type=independent".format(row['attributes'])]}
                            df_1 = pd.DataFrame(data = mydata)
                            rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)
                else:
                    sRNA_type.append("not_sRNA")
                    mydata = {'seqID': [self.genome_name], 'source': ['RendSeq'], 'feature': [row['feature']], 'start': [row['start']], 'end': [row['end']], 'score': [row['score']], 'strand': [row['strand']], 'phase': [row['phase']], 'attributes': ["{}".format(row['attributes'])]}
                    df_1 = pd.DataFrame(data = mydata)
                    rend_feature = pd.concat([rend_feature, df_1], ignore_index = True)

            rend_feature.to_csv(os.path.join(self.rend_folder, os.path.basename(file)), sep = '\t', index = False, header = None)


            if 'All_conditions' in file:
                statfile = self.analysis_folder + 'Statistics_All_conditions.txt'
                df = pd.read_csv(statfile, sep = '\t', header = 0)
                df['sRNA_type'] = sRNA_type
                df.to_csv(os.path.join(self.analysis_folder, 'Statistics_All_conditions.txt'), sep = '\t', index = False)
#
# # Figure 5 shows a screenshot of JBrowse that shows an example of what the annotation file produced by this
# # pipeline contains.
#
#
# ####################################################################
# ### Re-annotate operons as CDS, 5'UTR's, 3'UTRs, and intergenic RNAs.

    def categorize_operons(self):
        for file in glob.glob(self.rend_folder + '*.gff3'):
            df_rend = pd.read_csv(file, sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
            df_rend = df_rend[df_rend['feature'] == 'CDS']
            df_nic = self.df_Nicolas
            df_nic = df_nic[df_nic['feature'] == 'CDS']
            df_cat = pd.DataFrame()
            df_alternative = pd.DataFrame()
            for index, row in df_rend.iterrows():
                if 'Type=Alternative start site' in row['attributes']:
                    df_alternative = pd.concat([df_alternative, df_rend[df_rend.index == index]])
                else:
                    start = row['start']
                    end = row['end']
                    strand = row['strand']
                    condition = (df_nic['start'] >= start) & (df_nic['end'] <= end)
                    df_overlap = df_nic[condition][df_nic[condition]['strand'] == strand].sort_values(['start', 'end']).reset_index(drop = True)
                    if len(df_overlap) == 0:
                        pass
                    if len(df_overlap) > 0:
                        first_string = df_overlap['attributes'][0]
                        last_string = df_overlap['attributes'][max(df_overlap.index)]
                        first_list = first_string.split(';')
                        last_list = last_string.split(';')
                        for item in first_list:
                            if item.find('Name') == 0:
                                first_name = item
                        for item in last_list:
                            if item.find('Name') == 0:
                                last_name = item
                            try:
                                type(first_name)
                            except NameError:
                                first_name = 'ID{}..{}'.format(start, end)
                            try:
                                type(last_name)
                            except NameError:
                                last_name = 'ID{}..{}'.format(start, end)
                        if strand == '+':
                            UTR5 = pd.DataFrame({'seqID': [row['seqID']],'source': [row['source']],'feature': ["5'UTR"], 'start': [row['start']], 'end': [min(df_overlap['start'])],'score': [row['score']], 'strand': [row['strand']], 'phase': [row['phase']], 'attributes': ['{}'.format(first_name)]})
                            UTR3 = pd.DataFrame({'seqID': [row['seqID']],'source': [row['source']],'feature': ["3'UTR"], 'start': [max(df_overlap['end'])], 'end': [row['end']],'score': [row['score']], 'strand': [row['strand']], 'phase': [row['phase']], 'attributes': ['{}'.format(last_name)]})
                        else:
                            UTR3 = pd.DataFrame({'seqID': [row['seqID']],'source': [row['source']],'feature': ["3'UTR"], 'start': [row['start']], 'end': [min(df_overlap['start'])],'score': [row['score']], 'strand': [row['strand']], 'phase': [row['phase']], 'attributes': ['{}'.format(first_name)]})
                            UTR5 = pd.DataFrame({'seqID': [row['seqID']],'source': [row['source']],'feature': ["5'UTR"], 'start': [max(df_overlap['end'])], 'end': [row['end']],'score': [row['score']], 'strand': [row['strand']], 'phase': [row['phase']], 'attributes': ['{}'.format(last_name)]})
                    df_cat = pd.concat([df_cat, UTR5, UTR3])
                    if len(df_overlap) == 1:
                        CDS = pd.DataFrame({'seqID': [row['seqID']],'source': [row['source']],'feature': ['CDS'], 'start': [min(df_overlap['start'])], 'end': [max(df_overlap['end'])],'score': [row['score']], 'strand': [row['strand']], 'phase': [row['phase']], 'attributes': [df_overlap['attributes'][0]]})
                        df_cat = pd.concat([df_cat, CDS])
                    elif len(df_overlap) > 1:
                        for i in range(0, max(df_overlap.index)):
                            first_string = df_overlap['attributes'][i]
                            last_string = df_overlap['attributes'][i + 1]
                            first_list = first_string.split(';')
                            last_list = last_string.split(';')
                            for item in first_list:
                                if item.find('Name=') == 0:
                                    first_name = item
                            for item in last_list:
                                if item.find('Name=') == 0:
                                    last_name = item[5:]
                            if df_overlap['end'][i] < df_overlap['start'][i + 1]:
                                intergenic = pd.DataFrame({'seqID': [row['seqID']],'source': [row['source']],'feature': ["intergenic_UTR"], 'start': [df_overlap['end'][i]], 'end': [df_overlap['start'][i + 1]],'score': [row['score']], 'strand': [row['strand']], 'phase': [row['phase']], 'attributes': ['{}_{}'.format(first_name, last_name)]})
                            else:
                                pass
                            CDS = pd.DataFrame({'seqID': [row['seqID']],'source': [row['source']],'feature': ['CDS'], 'start': [df_overlap['start'][i]], 'end': [df_overlap['end'][i]],'score': [row['score']], 'strand': [row['strand']], 'phase': [row['phase']], 'attributes': [df_overlap['attributes'][i]]})
                            df_cat = pd.concat([df_cat, CDS])
                            try:
                                df_cat = pd.concat([df_cat, intergenic])
                            except:
                                pass
                        CDS = pd.DataFrame({'seqID': [row['seqID']],'source': [row['source']],'feature': ['CDS'], 'start': [df_overlap['start'][max(df_overlap.index)]], 'end': [df_overlap['end'][max(df_overlap.index)]],'score': [row['score']], 'strand': [row['strand']], 'phase': [row['phase']], 'attributes': [df_overlap['attributes'][max(df_overlap.index)]]})
                        df_cat = pd.concat([df_cat, CDS])

            df_rend = pd.read_csv(file, sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
            df_rend = df_rend[df_rend['feature'] != 'CDS']
            df_rend = pd.concat([df_rend, df_cat, df_alternative])
            df_rend = df_rend.drop_duplicates(subset = ['start', 'end', 'strand'], keep = 'first', ignore_index = True)
            df_rend = df_rend.sort_values(['start', 'end'], ascending = True)
            df_rend.to_csv(os.path.join(file + '_categorized_operons.gff3'), sep = '\t', index = False, header = None)

#
# ################################################################################
# # ### Find Nicolas annotations that are not covered in RendSeq
# # ### For coding annotations in Nicolas, both start and stop has to covered by RendSeq or else
# # ### we say that the annotation is missing
# # ### All other features has to be just partially covered

    def not_covered(self):

        rend_file = self.rend_folder + '/All_conditions.gff3'
        df_rend = pd.read_csv(rend_file, sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
        self.df_Nicolas = self.df_Nicolas.sort_values(['start', 'end'], axis = 0, ascending = True)
        self.df_Nicolas.reset_index(inplace = True, drop = True)
        df_new = pd.DataFrame()

        overlap = 0
        for index, row in self.df_Nicolas.iterrows():
            N_strand = row['strand']
            N_start = row['start']
            N_end = row['end']
            ### If there is an overlap then one of these conditions will provide a dataframe in df_rend with a length above 0
            ### This method can be discussed as a RendSeq annotation might overlap with a Nicolas et al. annotation even though they are not the same genes
            ### e.g. yabG gene
            conditions = (N_start <= df_rend['start']) & (N_end > df_rend['start']) | (N_start >= df_rend['start']) & (N_start < df_rend['end']) | (N_start == df_rend['start']) & (N_end == df_rend['end'])
            if len(df_rend[conditions][df_rend[conditions]['strand'] == row['strand']]) > 0:
                overlap += 1
            else:
                ### Make a dataframe with annotations not covered by Rend-seq
                df_new = pd.concat([df_new, self.df_Nicolas[self.df_Nicolas.index == index]], ignore_index = True)

        ### Save annotations not covered by Rend-seq in file
        df_new.to_csv(os.path.join(self.analysis_folder, 'annotations_not_in_RendSeq.gff3'), sep = '\t', index = False, header = None)
#
# ################################################################################
# ########### Supply missing annotations with Nicolas et al. #####################
# ### For genes not covered by RendSeq supply with Nicolas et al.
# ### CDS need to be to found in Rend-seq file at exact same positions, not just a partial overlap
# ### All other annotations need just a partial overlap. If no overlap we transfer gene from Nicolas et al.

    def supply_cat_operon(self):
        for file in glob.glob(self.rend_folder + 'All_conditions.gff3_categorized_operons.gff3'):

            df_rend = pd.read_csv(file, sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
            df_rend = df_rend[df_rend['feature'] == 'CDS']
            df_nic = self.df_Nicolas
            df_nic = df_nic[df_nic['feature'] == 'CDS']

            Rend_seq_start = df_rend['start'][df_rend['feature'] == 'CDS'].tolist()
            Rend_seq_end = df_rend['end'][df_rend['feature'] == 'CDS'].tolist()
            Rend_seq_strand = df_rend['strand'][df_rend['feature'] == 'CDS'].tolist()
            Rend_seq_tupples = []
            for start, end, strand in zip(Rend_seq_start, Rend_seq_end, Rend_seq_strand):
                Rend_seq_tupples.append((start, end, strand))

            nicolas_start = df_nic['start'][df_nic['feature'] == 'CDS'].tolist()
            nicolas_end = df_nic['end'][df_nic['feature'] == 'CDS'].tolist()
            nicolas_strand = df_nic['strand'][df_nic['feature'] == 'CDS'].tolist()
            nicolas_index = df_nic.index[df_nic['feature'] == 'CDS'].tolist()
            nicolas_tupples = []
            for start, end, strand, index in zip(nicolas_start, nicolas_end, nicolas_strand, nicolas_index):
                nicolas_tupples.append((start, end, strand, index))

            missing_list = []
            for item in nicolas_tupples:
                tupple = (item[0], item[1], item[2])
                if tupple not in Rend_seq_tupples:
                     missing_list.append(item[3])

            df_rend = pd.read_csv(file, sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
            for item in missing_list:
                df_rend = pd.concat([df_rend, df_nic[df_nic.index == item]], ignore_index = True)

            ### now add all other features that are not CDS and that don't overlap with any Rend-seq annotations
            df_nic_no_overlap = self.analysis_folder + 'annotations_not_in_RendSeq.gff3'
            df_nic_no_overlap = pd.read_csv(df_nic_no_overlap, sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

            df_rend = pd.concat([df_rend, df_nic_no_overlap[df_nic_no_overlap['feature'] != 'CDS']], ignore_index = True)
            df_rend = df_rend.sort_values(['start', 'end'], axis = 'index', ascending = True)

            df_rend.to_csv(os.path.join(file + '_supplied_with_{}.gff3'.format(os.path.basename(self.old_annotation_name))), sep = '\t', index = False, header = None)

    def supply_raw(self):
        df_rend = pd.read_csv(self.rend_folder + 'All_conditions.gff3', sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
        df_operon = pd.read_csv(self.rend_folder + 'All_conditions.gff3_categorized_operons.gff3' + '_supplied_with_{}.gff3'.format(os.path.basename(self.old_annotation_name)), sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
        df_operon = df_operon[(df_operon['source'] != 'RendSeq') & (df_operon['feature'] == 'CDS')]
        df_rend = pd.concat([df_rend, df_operon])
        df_nic_no_overlap = self.analysis_folder + 'annotations_not_in_RendSeq.gff3'
        df_nic_no_overlap = pd.read_csv(df_nic_no_overlap, sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

        df_rend = pd.concat([df_rend, df_nic_no_overlap[df_nic_no_overlap['feature'] != 'CDS']], ignore_index = True)
        df_rend = df_rend.sort_values(['start', 'end'], axis = 'index', ascending = True)

        df_rend.to_csv(os.path.join(self.rend_folder, 'All_conditions.gff3_' + 'supplied_with_{}.gff3'.format(os.path.basename(self.old_annotation_name))), sep = '\t', index = False, header = None)


def delete_1_bp(file):
    df = pd.read_csv(file, sep = '\t', names = ['seqID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    for index, row in df.iterrows():
        if abs(row['start'] - row['end']) <= 1:
            df = df[df.index != index]
    df.to_csv(file, sep = '\t', index = False, header = None)
