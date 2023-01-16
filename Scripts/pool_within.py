import os
import pandas as pd
import glob
import statistics as stat

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


    # # # ###########################################################################################
    # # # ########## Pool annotations where start and stop sites are within 5 bp #####################
    # # ### The pooled dataframe with annotations from all conditions is being compared to itself
    # # ### For any annotation we check if it overlaps with any other annotations and if so
    # # ### if the distance between start and stop site is within 5 bp.
    # # # ### Annotations that have the same start and stop (+/- 5 bp) will be pooled into óne annotation
    # # # ### Whatever start position is the most frequent occuring will be set as canonical and in attributes it will be
    # # # ### given the distance from the actual position of any condition to the canonical.
    # # ### We then count how many annotations where this is true and divide by the total number of annotations to
    # # ### provide a score of the fraction of conditions where the annotation is present.


def pool(file, distance, wig_folder, output_folder, genome_name):
    df_rend = pd.read_csv(file, sep = '\t', header = 0)
    df_rend.reset_index(inplace = True, drop = True)
    df_new = pd.DataFrame({'seqID': [], 'source': [], 'feature': [], 'start': [], 'end': [], 'score': [], 'strand': [], 'phase': [], 'attributes': []})
    rend_startlist = df_rend['start'].tolist()
    rend_endlist = df_rend['end'].tolist()
    rend_strandlist = df_rend['strand'].tolist()
    rend_attributelist = df_rend['attributes'].tolist()
    rend_phaselist = df_rend['phase'].tolist()

    print('Pool annotations within {} bp of each other'.format(distance))

    for i in range(0,len(rend_startlist)):
        StartR = rend_startlist[i]
        EndR = rend_endlist[i]
        StrandR = rend_strandlist[i]
        df_copy = df_rend[df_rend['strand'] == StrandR].copy()
        ### conditions to check if the annotation number i on the list overlaps with any other annotions
        conditions = (StartR < df_copy['end']) & (df_copy['end'] <= EndR) | (StartR <= df_copy['start']) & (EndR >= df_copy['end']) | (StartR <= df_copy['start']) & (df_copy['start'] < EndR) | (df_copy['start'] < StartR) & (df_copy['end'] > EndR)
        ### if it overlaps with itself only then just add it to the pooled dataframe
        if len(df_copy[conditions]) == 1:
            mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [rend_startlist[i]], 'end': [rend_endlist[i]], 'score': [1], 'strand': [rend_strandlist[i]], 'phase': [rend_phaselist[i]], 'attributes': ['condition={}'.format(rend_attributelist[i])]}
            df_1 = pd.DataFrame(data = mydata)
            df_new = pd.concat([df_new, df_1], ignore_index = True)
        else:
            startlist = []
            endlist = []
            phaselist = []
            df_copy = df_copy[conditions]
            copy_strand = df_copy['strand'].tolist()
            copy_start = df_copy['start'].tolist()
            copy_end = df_copy['end'].tolist()
            copy_phase = df_copy['phase'].tolist()
            for n in range(0, len(df_copy)):
                StrandL = copy_strand[n]
                StartL = copy_start[n]
                EndL = copy_end[n]
                if abs(StartR - StartL) <= distance and abs(EndR - EndL) <= distance:
                    startlist.append(StartL)
                    endlist.append(EndL)
                    phaselist.append(copy_phase[n])
                else:
                    pass
            if len(startlist) == 1:
                thephase = rend_phaselist[i]
            elif len(startlist) > 1 and '.' in phaselist:
                thephase = '.'
            else:
                thephase = phaselist[0]
            ### if the only annotation that lies within distance bp is itself
            if len(startlist) == 1:
                mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [rend_startlist[i]], 'end': [rend_endlist[i]], 'score': [1], 'strand': [rend_strandlist[i]], 'phase': [rend_phaselist[i]], 'attributes': ['condition={}'.format(rend_attributelist[i])]}
                df_1 = pd.DataFrame(data = mydata)
                df_new = pd.concat([df_new, df_1], ignore_index = True)
            else:
                ### we will use the most frequent start and stop site as cannonical in the pooled dataframe (df_new), so if we subset the pooled dataframe using those conditions
                ### and retrieve and empty dataframe, then this current annotation is the first to be added.
                conditions2 = (df_new['start'] == stat.mode(startlist)) & (df_new['end'] == stat.mode(endlist)) & (df_new['strand'] == StrandR)
                if len(df_new[conditions2]) == 0:
                    if StrandR == '+':
                        startdist = StartR - stat.mode(startlist)
                        enddist = EndR - stat.mode(endlist)
                    else:
                        startdist = stat.mode(startlist) - StartR
                        enddist = stat.mode(endlist) - EndR
                    ### make a row in the pooled dataframe using the most frequent occuring start and stop sites and add to attributes the condition of the
                    ### annotation in question (its attribute). Further add to attributes the distance from the current annotation the one in the pooled dataframe
                    ### distance is given by setting the cannonical site as 0 and upstream from that is + and downstream is -
                    mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [stat.mode(startlist)], 'end': [stat.mode(endlist)], 'score': [1], 'strand': [StrandR], 'phase': [thephase], 'attributes': ['{}=start:{}, stop:{}'.format(rend_attributelist[i], startdist, enddist)]}
                    df_1 = pd.DataFrame(data = mydata)
                    df_new = pd.concat([df_new, df_1], ignore_index = True)
                ### being here means that the current annotation is the not the first to be added to the pooled dataframe and hence we will add the current attribute to the
                ### attribute of the pooled annotation. We also add a score of 1 for every time we add another annotation.
                else:
                    if StrandR == '+':
                        startdist = StartR - stat.mode(startlist)
                        enddist = EndR - stat.mode(endlist)
                    else:
                        startdist = stat.mode(startlist) - StartR
                        enddist = stat.mode(endlist) - EndR
                    mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [stat.mode(startlist)], 'end': [stat.mode(endlist)], 'score': [1 + df_new[conditions2]['score'].iloc[0]], 'strand': [StrandR], 'phase': [thephase], 'attributes': ['{};{}=start:{}, stop:{}'.format(df_new[conditions2]['attributes'].iloc[0], rend_attributelist[i], startdist, enddist)]}
                    df_1 = pd.DataFrame(data = mydata)
                    df_new = df_new.drop(df_new[conditions2].index.values[0], axis = 0)
                    df_new = pd.concat([df_new, df_1], ignore_index = True)



    thescore = df_new['score']
    newscore = []
    ### by dividing the score of each annotation in the pooled dataframe by the total number of conditions
    ### we get the fraction of samples of which the given annotation is expressed.
    for item in thescore:
        newscore.append(item / (len(glob.glob(wig_folder + '/*.wig')) / 4))

    df_new['score'] = newscore
    df_new.to_csv(os.path.join(output_folder, 'All_conditions.gff3'), sep = '\t', index = False)
