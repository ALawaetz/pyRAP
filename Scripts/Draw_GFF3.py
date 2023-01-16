import pandas as pd
import os
from multiprocessing import Process, Pool
import folders
import sys
import glob

input_folder = sys.argv[1]
output_folder = sys.argv[2]
operon_file = sys.argv[3]
genome_name = sys.argv[4]
primary_operons = sys.argv[5]




## We use strict parameters to determine starts and stops of Rend-seq operons (read density 0.8) but some genes
## are expressed at so low levels that the read density is below 0.8 and those genes will therefore not be assigned
## an operon with above method. In such cases we assign to BSGatlas operons.
## The BSG atlas (https://rth.dk/resources/bsgatlas/index.php) provide an annotation file
## of operons in Bacillus subtilis.

################################################################################
################################################################################
# #################### draw GFF3 ##########################
# ## We now have all the processed files necessary to make an annotation file in GFF3 format.

### With the two operon files we assign each peak in the processed wig files to an operon
### first if possible (overlap), we assing to Rend-seq operons. For remaining peaks that have not been assigned,
### we assign if possible (overlap) to BSG operons. Still after that, some peaks fall outside of operons, and they will be
### named by the nearest BSG atlas operon and added 'a_half' to their name
### By such we can afterwards connect all start sites to all stop sites within the same operon.

# # If a start doesn't overlap with an operon, then see if the downstream stop overlaps with one and then
# # assign that operon to both peaks.
# # If a stop doesn't overlap with an operon, then see if the upstream start overlaps with one and then
# # assign that operon to both peaks.
# # If none of them overlaps with an operon then name them by the upstream operon + '_and_a_half', unless there is an operon between the two peaks,
# # in which case the peaks will be assigned to that operon.
# # When drawing annotations, draw annotations from each start site to each downstream stop site that belongs to the same operon and
# # where the read density between peaks is above 0.2 (arbitrary number)

def YourCode(file):
    if '.csv' in file and 'GFF3' not in file:
        df_rend = pd.read_csv(file, sep = '\t', header = 0)
        try:
            thestrand = max(df_rend['strand'])
        except ValueError:
            print('''No annotations made for {}
            Try analysing a larger region or lowering your cutoff values.
            For help type "python3 pipeline.py --help"'''.format(os.path.basename(file)))
        df_rend = df_rend[df_rend['peak'] != '.']
        df_rend.reset_index(drop = True, inplace = True)
        rend_location = df_rend['location'].tolist()
        rend_peak = df_rend['peak'].tolist()
        del df_rend

        if primary_operons == 'None':
            operonfile = glob.glob(folders.operon_folder + '*_FwdandRev.gff3')
            df_operon = pd.read_csv(operonfile[0], sep = '\t', usecols = ['start', 'end', 'strand', 'attributes'])
            df_operon = df_operon.sort_values(['strand', 'start'], ascending = True, axis = 0)
        else:
            df_operon = pd.read_csv(primary_operons, sep = '\t', names = ['seqID', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
            df_operon = df_operon.sort_values(['strand', 'start'], ascending = True, axis = 0)
        ### make an empty list to add operon information

        attributes = []

        ### Divide operon file into fwd and rev operons so there is less to iterate over.

        if len(rend_peak) > 0 and primary_operons == 'None':
            if thestrand == 0:
                df_operon = df_operon[df_operon['strand'] == 0]
            else:
                df_operon = df_operon[df_operon['strand'] == 1]
                df_operon = df_operon.sort_values(['end'], ascending = False, axis = 0)
        elif len(rend_peak) > 0 and primary_operons != 'None':
            if thestrand == 0:
                df_operon = df_operon[df_operon['strand'] == '+']
            else:
                df_operon = df_operon[df_operon['strand'] == '-']
                df_operon = df_operon.sort_values(['end'], ascending = False, axis = 0)
        else:
            pass
        df_operon.reset_index(drop = True)
        operon_start = df_operon['start'].tolist()
        operon_end = df_operon['end'].tolist()
        operon_attributes = df_operon['attributes'].tolist()
        del df_operon

    ## Peaks not assigned to Rend_seq operons will be assigned 'Cuckoo' which will be attempted to assign to BSG operons followingly
        for i in range(0, len(rend_peak)):
            x = 1
            ### If the peak is a START then iterate bottom up ([::-1])
            ### else iterate top down
            if rend_peak[i] == 'START':
                for z in range(0, len(operon_start)):
                    ### if there is an overlap with any operons then add that operon to the list
                    if rend_location[i] >= operon_start[z] and rend_location[i] <= operon_end[z]:
                        attributes.append(operon_attributes[z])
                        break
                    else:
                        ### if the start doesn't overlap with any operons
                        ### and all rows in df_operon has been itereated
                        ### see if the downstream stop overlaps with any
                        if x == len(operon_start):
                            x = 1
                            n = 0
                            ### find the next downstream stop and use its location to iterate through df_operon
                            while rend_peak[i + n] == 'START':
                                n += 1
                                if i + n > len(rend_peak) - 1:
                                    attributes.append('Cuckoo')
                                    break
                                else:
                                    pass
                                if rend_peak[i + n] == 'STOP':
                                    for t in range(0, len(operon_start)):
                                        ### if this peak overlaps with an operon or if there is an operon inbetween this peak and the upstream start-peak
                                        ### Then append that operon to the list.
                                        if t != 0:
                                            if thestrand == 0 and rend_location[i + n] >= operon_start[t] and rend_location[i + n] <= operon_end[t] and operon_end[t - 1] < rend_location[i] or thestrand == 0 and rend_location[i] <= operon_start[t] and rend_location[i + n] >= operon_end[t] and operon_end[t - 1] < rend_location[i] or thestrand == 1 and rend_location[i + n] >= operon_start[t] and rend_location[i + n] <= operon_end[t] and operon_end[t - 1] > rend_location[i] or thestrand == 1 and rend_location[i] >= operon_end[t] and rend_location[i + n] <= operon_start[t] and operon_end[t - 1] > rend_location[i]:
                                                attributes.append(operon_attributes[t])
                                                break
                                            else:
                                                if x == len(operon_start):
                                                    ### being here means that the peak is between operons.
                                                    ### so we will name after the upstream operon + '_and_a_half'
                                                    ### that will be the last operon in the operon dataframe when we include only the rows
                                                    ### where 'end' is less than the position of the peak
                                                    attributes.append('Cuckoo')
                                                    break
                                                else:
                                                    x += 1
                                        else:
                                            if thestrand == 0 and rend_location[i + n] >= operon_start[t] and rend_location[i + n] <= operon_end[t] or thestrand == 0 and rend_location[i] <= operon_start[t] and rend_location[i + n] >= operon_end[t] or thestrand == 1 and rend_location[i + n] >= operon_start[t] and rend_location[i + n] <= operon_end[t] or thestrand == 1 and rend_location[i] >= operon_end[t] and rend_location[i + n] <= operon_start[t]:
                                                attributes.append(operon_attributes[t])
                                                break
                                            else:
                                                if x == len(operon_start):
                                                    ### being here means that the peak is between operons.
                                                    ### so we will name after the upstream operon + '_and_a_half'
                                                    ### that will be the last operon in the operon dataframe when we include only the rows
                                                    ### where 'end' is less than the position of the peak
                                                    attributes.append('Cuckoo')
                                                    break
                                                else:
                                                    x += 1
                                    break
                                else:
                                    pass

                        else:
                            ### make a count so we now when we have iterated on all rows
                            x += 1
            else:
                ### this is all the same as above except if a stop doesn't overlap with any operons
                ### we see if the upstream start does. If not, add 'Cuckoo' to the list.
                for z in range(0, len(operon_start)):
                    if rend_location[i] >= operon_start[z] and rend_location[i] <= operon_end[z]:
                        attributes.append(operon_attributes[z])
                        break
                    else:
                        if x == len(operon_start):
                            x = 1
                            n = 0
                            while rend_peak[i - n] == 'STOP':
                                n += 1
                                if i - n < 0:
                                    attributes.append('Cuckoo')
                                    break
                                else:
                                    pass
                                if rend_peak[i - n] == 'START':
                                    for t in range(0, len(operon_start)):
                                        ### if this peak overlaps with an operon or if there is an operon inbetween this peak and the downstream stop-peak
                                        ### Then append that operon to the list.
                                        if t + 1 != len(operon_start):
                                            if thestrand == 0 and rend_location[i - n] >= operon_start[t] and rend_location[i - n] <= operon_end[t] and operon_start[t + 1] > rend_location[i] or thestrand == 0 and rend_location[i] >= operon_end[t] and rend_location[i - n] <= operon_start[t] and operon_start[t + 1] > rend_location[i] or thestrand == 1 and rend_location[i - n] >= operon_start[t] and rend_location[i - n] <= operon_end[t] and operon_start[t + 1] < rend_location[i] or thestrand == 1 and rend_location[i - n] >= operon_end[t] and rend_location[i] <= operon_start[t] and operon_start[t + 1] < rend_location[i]:
                                                attributes.append(operon_attributes[t])
                                                break
                                            else:
                                                if x == len(operon_start):
                                                    attributes.append('Cuckoo')
                                                    break
                                                else:
                                                    x += 1
                                        else:
                                            if thestrand == 0 and rend_location[i - n] >= operon_start[t] and rend_location[i - n] <= operon_end[t] or thestrand == 0 and rend_location[i] >= operon_end[t] and rend_location[i - n] <= operon_start[t] or thestrand == 1 and rend_location[i - n] >= operon_start[t] and rend_location[i - n] <= operon_end[t] or thestrand == 1 and rend_location[i - n] >= operon_end[t] and rend_location[i] <= operon_start[t]:
                                                attributes.append(operon_attributes[t])
                                                break
                                            else:
                                                if x == len(operon_start):
                                                    attributes.append('Cuckoo')
                                                    break
                                                else:
                                                    x += 1
                                    break
                                else:
                                    pass


                        else:
                            x += 1

        ### add the operon data to the original dataframe
        ###
        ### read the file again, but use all columns and then add the operon list
        df_rend = pd.read_csv(file, sep = '\t', header = 0)
        df_rend = df_rend[df_rend['peak'] != '.']
        if len(attributes) != 0:
            df_rend['operon'] = attributes
        else:
            df_rend['operon'] = ['.'] * len(df_rend)

        ### read the file once more, but now include all rows
        df_rend_full = pd.read_csv(file, sep = '\t', header = 0)

        ### concat the two dataframes and then drop the duplicates and keep the last added
        df_rend = pd.concat([df_rend_full, df_rend])
        df_rend = df_rend.fillna('.')
        df_rend = df_rend.drop_duplicates(subset = ['location', 'end'], keep = 'last')

        if thestrand == 0:
            df_rend = df_rend.sort_values(['location', 'end'], axis = 0, ascending = True)
        else:
            df_rend = df_rend.sort_values(['location', 'end'], axis = 0, ascending = False)

        df_rend.to_csv(file, sep = '\t', index = False)

        del df_rend_full
    #
    ### now, for all the 'Cuckoo's assign them to a BSG-operon
    ### Method is bacically the same as above for Rend-seq operons

        if operon_file == 'None':
            del df_rend
        else:
            thestrand = df_rend['strand'][0]
            df_rend = df_rend[df_rend['peak'] != '.']
            df_rend.reset_index(drop = True, inplace = True)
            rend_location = df_rend['location'].tolist()
            rend_peak = df_rend['peak'].tolist()
            rend_operon = df_rend['operon'].tolist()
            del df_rend
            df_operon = pd.read_csv(operon_file, sep = '\t', usecols = ['start', 'end', 'strand', 'attributes'])
            df_operon = df_operon.sort_values(['strand', 'start'], ascending = True, axis = 0)

            ### make an empty list to add operon information
            attributes = []

            ### Remove excess rows so there is less to iterate over.
            if len(rend_peak) > 0:
                if thestrand == 0:
                    df_operon = df_operon[df_operon['strand'] == '+']
                else:
                    df_operon = df_operon[df_operon['strand'] == '-']
                    ### Some operons overlap
                    ### if a start site overlaps with two operons it should be assigned to the second operon on the list which have
                    ### been ordered by start position in an ascending fascion. Or if we go bottom-up from the list, it should be
                    ### the first operon on the list.
                    ### if a stop site overlaps with two operons it should be assigned to the first operon on the list which have
                    ### been ordered by start position in an ascending fascion.
                    df_operon = df_operon.sort_values(['end'], ascending = False, axis = 0)
            else:
                pass

            df_operon.reset_index(drop = True, inplace = True)
            operon_start = df_operon['start'].tolist()
            operon_end = df_operon['end'].tolist()
            operon_attributes = df_operon['attributes'].tolist()
            del df_operon

            for i in range(0, len(rend_peak)):
                x = 1
                ### If the peak is a START then iterate bottom up ([::-1])
                ### else iterate top down
                if rend_operon[i] != 'Cuckoo':
                    attributes.append(rend_operon[i])
                elif rend_peak[i] == 'START': #and rend_peak[i] == 'Cuckoo':
                    operon_start = operon_start[::-1]
                    operon_end = operon_end[::-1]
                    operon_attributes = operon_attributes[::-1]
                    for z in range(0, len(operon_start)):
                        ### if there is an overlap with any operons then add that operon to the list
                        if rend_location[i] >= operon_start[z] and rend_location[i] <= operon_end[z]:
                            attributes.append(operon_attributes[z])
                            break
                        else:
                            ### if the start doesn't overlap with any operons
                            ### and all rows in df_operon has been itereated
                            ### see if the downstream stop overlaps with any
                            if x == len(operon_end):
                                x = 1
                                n = 0
                                ### find the next downstream stop and use its location to iterate through df_operon
                                while rend_peak[i + n] == 'START':
                                    n += 1
                                    if i + n > len(rend_peak) - 1 or rend_peak[i + n] == 'STOP':
                                        for t in range(0, len(operon_start)):
                                            ### if this peak overlaps with an operon or if there is an operon inbetween this peak and the upstream start-peak
                                            ### Then append that operon to the list.
                                            if i + n <= len(rend_peak) - 1 and rend_location[i + n] >= operon_start[t] and rend_location[i + n] <= operon_end[t] or i + n <= len(rend_peak) - 1 and thestrand == 0 and rend_location[i] <= operon_start[t] and rend_location[i + n] >= operon_end[t] or i + n <= len(rend_peak) - 1 and thestrand == 1 and rend_location[i] >= operon_end[t] and rend_location[i + n] <= operon_start[t]:
                                                attributes.append(operon_attributes[t])
                                                break
                                            else:
                                                ### if we have iterated through all rows and there is no overlap
                                                ### then add 'Cuckoo' to the list
                                                if x == len(operon_start):
                                                    ### being here means that the peak is between operons.
                                                    ### so we will name after the upstream operon (if any) + '_and_a_half'
                                                    ### that will be the last operon in the operon dataframe when we include only the rows
                                                    ### where 'end' is less than the position of the peak
                                                    if thestrand == 0:
                                                        if len([item for item in operon_end if item < rend_location[i]]) != 0:
                                                            ### the lenth of the attributes list be the same length as the operon list which can be indiced as the line above, and then we take
                                                            ### the last item on that list.
                                                            attributes.append(operon_attributes[len(operon_end) - len([item for item in operon_end if item < rend_location[i]]):][0] + '_an_a_half')
                                                            #attributes.append(df_operon[df_operon['end'] < rowRend['location']][::-1].reset_index()['attributes'][0] + '_an_a_half')
                                                            break
                                                        ### if there is no upstream operon we will name it by its downstream operon + '_minus_a_half'
                                                        ### which will be the first row in the operon dataframe when we include only the rows
                                                        ### where start is more than the position of the peak
                                                        else:
                                                            ### unless it is on the reverse strand in which case it is the first item on the list
                                                            attributes.append(operon_attributes[-1] + '_minus_a_half')
                                                            break
                                                    ### with reverse strand files we have to do it opposite
                                                    else:
                                                        if len([item for item in operon_end if item > rend_location[i]]) != 0:
                                                            attributes.append(operon_attributes[len(operon_end) - len([item for item in operon_end if item > rend_location[i]]):][0] + '_an_a_half')
                                                            break
                                                        else:
                                                            attributes.append(operon_attributes[-1] + '_minus_a_half')
                                                            break
                                                else:
                                                    x += 1
                                        break
                                    else:
                                        pass
                            else:
                                ### make a count so we know when we have iterated on all rows
                                x += 1
                    operon_start = operon_start[::-1]
                    operon_end = operon_end[::-1]
                    operon_attributes = operon_attributes[::-1]
                elif rend_peak[i] == 'STOP':# and rowRend['operon'] == 'Cuckoo':
                    ### this is all the same as above except if a stop doesn't overlap with any operons
                    ### we see if the upstream start does. If not, add 'Cuckoo' to the list.
                    operon_start = operon_start[::-1]
                    operon_end = operon_end[::-1]
                    operon_attributes = operon_attributes[::-1]
                    for z in range(0, len(operon_start)):
                        if rend_location[i] >= operon_start[z] and rend_location[i] <= operon_end[z]:
                            attributes.append(operon_attributes[z])
                            break
                        else:
                            if x == len(operon_end):
                                x = 1
                                n = 0
                                while rend_peak[i - n] == 'STOP':
                                    n += 1
                                    if i - n < 0 or rend_peak[i - n] == 'START':
                                        for t in range(0, len(operon_start)):
                                            ### if this peak overlaps with an operon or if there is an operon inbetween this peak and the downstream stop-peak
                                            ### Then append that operon to the list.
                                            if i - n >= 0 and rend_location[i - n] >= operon_start[t] and rend_location[i - n] <= operon_end[t] or i - n >= 0 and thestrand == 0 and rend_location[i] >= operon_end[t] and rend_location[i - n] <= operon_start[t] or i - n >= 0 and thestrand == 1 and rend_location[i] <= operon_start[t] and rend_location[i - n] >= operon_end[t]:
                                                attributes.append(operon_attributes[t])
                                                break
                                            else:
                                                if x == len(operon_end):
                                                    if thestrand == 0:
                                                        if len([item for item in operon_end if item < rend_location[i]]) != 0:
                                                            ### the lenth of the attributes list be the same length as the operon list which can be indiced as the line above, and then we take
                                                            ### the last item on that list.
                                                            attributes.append(operon_attributes[len(operon_end) - len([item for item in operon_end if item < rend_location[i]]):][0] + '_an_a_half')
                                                            #attributes.append(df_operon[df_operon['end'] < rowRend['location']][::-1].reset_index()['attributes'][0] + '_an_a_half')
                                                            break
                                                        ### if there is no upstream operon we will name it by its downstream operon + '_minus_a_half'
                                                        ### which will be the first row in the operon dataframe when we include only the rows
                                                        ### where start is more than the position of the peak
                                                        else:
                                                            ### unless it is on the reverse strand in which case it is the first item on the list
                                                            attributes.append(operon_attributes[-1] + '_minus_a_half')
                                                            break
                                                    ### with reverse strand files we have to do it opposite
                                                    else:
                                                        if len([item for item in operon_end if item > rend_location[i]]) != 0:
                                                            attributes.append(operon_attributes[len(operon_end) - len([item for item in operon_end if item > rend_location[i]]):][0] + '_an_a_half')
                                                            break
                                                        else:
                                                            attributes.append(operon_attributes[-1] + '_minus_a_half')
                                                            break
                                                else:
                                                    x += 1
                                        break
                                    else:
                                        pass


                            else:
                                x += 1
                else:
                    pass

            ### add the operon data to the original dataframe
            ###
            ### read the file again, but use all columns and then add the operon list
            df_rend = pd.read_csv(file, sep = '\t', header = 0)
            df_rend = df_rend[df_rend['peak'] != '.']

            df_rend['operon'] = attributes

            ### read the file once more, but now include all rows
            df_rend_full = pd.read_csv(file, sep = '\t', header = 0)

            ### concat the two dataframes and then drop the duplicates and keep the last added
            df_rend = pd.concat([df_rend_full, df_rend])
            df_rend = df_rend.fillna('.')
            df_rend = df_rend.drop_duplicates(subset = ['location', 'end'], keep = 'last')

            if thestrand == 0:
                df_rend = df_rend.sort_values(['location', 'end'], axis = 0, ascending = True)
            else:
                df_rend = df_rend.sort_values(['location', 'end'], axis = 0, ascending = False)

            df_rend.to_csv(file, sep = '\t', index = False)

            del df_rend, df_rend_full

        ### Draw GFF3 annotation
        name = file
        myheader = {'seqID': [], 'source': [], 'feature': [], 'start': [], 'end': [], 'score': [], 'strand': [], 'phase': [], 'attributes': []}
        GFF3 = pd.DataFrame(data = myheader)

        ### Load dataframe in df_full and omit all rows not containing peaks and save as df
        ### we will iterate through each row in df
        ### and to asses the read density between peaks we cross reference to df_full
        df_full = pd.read_csv(file, sep = '\t', usecols = ['location', 'value', 'strand', 'peak', 'operon'])
        thestrand = df_full['strand'][0]
        values = df_full['value'].tolist()
        full_location = df_full['location'].tolist()

        df = df_full[df_full['peak'] != '.']
        #df = df.reset_index(drop = True)

        df_location = df['location'].tolist()
        df_peak = df['peak'].tolist()
        df_operon = df['operon'].tolist()
        del df, df_full

        for i in range(0, len(df_peak)):
            if df_peak[i] == 'START':
                read_dens = 1
                moving_dens = 1
                n = 0
                while i < len(df_peak) - 1 and read_dens >= 0.2 and moving_dens >= 0.1 and df_operon[i] == df_operon[i + n + 1]:

                    ### also make sure that the first 25 bp downstream a start and upstream a stop has density above 0.2
                    index_start = full_location.index(df_location[i + n])
                    index_end = full_location.index(df_location[i + n + 1])
                    condition_values = values[index_start: index_end + 1]
                    if abs(df_location[i + n] - df_location[i+ n + 1]) > 100:
                        for c in range(0, len(condition_values), 100):
                            for_condition_values = condition_values[c: c + 100]
                            non_zero_values = [z for z in for_condition_values if z != 0]
                            moving_dens = len(non_zero_values)/len(for_condition_values)
                            if moving_dens < 0.1:
                                break
                            else:
                                pass
                    else:
                        non_zero_values = [z for z in condition_values if z != 0]
                        read_dens = len(non_zero_values)/len(condition_values)

                    if read_dens >= 0.2 and moving_dens >= 0.1 and df_peak[i + n + 1] == 'STOP':
                        if thestrand == 0:
                            mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df_location[i]], 'end': [df_location[i + n + 1]], 'score': [1000000], 'strand': ['+'], 'phase': ['.'], 'attributes': [name]}
                        else:
                            mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df_location[i + n + 1]], 'end': [df_location[i]], 'score': [1000000], 'strand': ['-'], 'phase': ['.'], 'attributes': [name]}
                        newrow = pd.DataFrame(data = mydata)
                        GFF3 = pd.concat([GFF3, newrow])
                    else:
                        pass
                    n += 1
                    if i + n == len(df_peak) - 1:
                        break
                    else:
                        pass
            else:
                pass

        ### now some peaks will still not have been connected in transcript because they been assigned to different operons
        ### Those such peaks will be connected to the nearest corresponding start or stop if the read density is above 0.2
        ### and only if one of the peaks is assigned to an operon with the ending 'a_half' (peaks that falled outside both
        ### rend-seq and bsg operons). This means that peaks from different operons will still not be connected.

        GFF3_remain = pd.DataFrame(data = myheader)

        ### Load dataframe in df_full and omit all rows not containing peaks and save as df
        ### we will iterate through each row in df
        ### and to asses the read density between peaks we cross reference to df_full

        startlist = GFF3['start'][GFF3['strand'] == '+'].tolist()
        startlist = startlist + GFF3['end'][GFF3['strand'] == '-'].tolist()
        endlist = GFF3['end'][GFF3['strand'] == '+'].tolist()
        endlist = endlist + GFF3['start'][GFF3['strand'] == '-'].tolist()
        for i in range(0, len(df_peak)):
            if df_peak[i] == 'START' and df_location[i] not in startlist:
                read_dens = 1
                moving_dens = 1
                n = 0
                ### make a empty list to append operons. When connecting starts to stops
                ### We can move from operon X to operon X__and_a_half but not to operon X+1
                move_operon = []
                move_operon.append(df_operon[i])
                while i < len(df_peak) - 1 and read_dens >= 0.2 and moving_dens >= 0.1:

                    index_start = full_location.index(df_location[i + n])
                    index_end = full_location.index(df_location[i + n + 1])
                    condition_values = values[index_start: index_end + 1]
                    if abs(df_location[i + n] - df_location[i+ n + 1]) > 200:
                        for c in range(0, len(condition_values), 200):
                            for_condition_values = condition_values[c: c + 200]
                            non_zero_values = [z for z in for_condition_values if z != 0]
                            moving_dens = len(non_zero_values)/len(for_condition_values)
                            if moving_dens < 0.1:
                                break
                            else:
                                pass
                    else:
                        non_zero_values = [z for z in condition_values if z != 0]
                        read_dens = len(non_zero_values)/len(condition_values)

                    if read_dens >= 0.2 and moving_dens >= 0.1 and df_peak[i + n + 1] == 'STOP':
                        if '_a_half' in df_operon[i] or '_a_half' in df_operon[i + n + 1]:
                            if thestrand == 0:
                                mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df_location[i]], 'end': [df_location[i + n + 1]], 'score': [1000000], 'strand': ['+'], 'phase': ['.'], 'attributes': [name]}
                            else:
                                mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df_location[i + n + 1]], 'end': [df_location[i]], 'score': [1000000], 'strand': ['-'], 'phase': ['.'], 'attributes': [name]}
                            newrow = pd.DataFrame(data = mydata)
                            GFF3_remain = pd.concat([GFF3_remain, newrow])
                            break
                        else:
                            move_operon.append(df_operon[i + n + 1])
                            if len(list(set(move_operon))) >= 2:
                                break
                            else:
                                pass
                    else:
                        pass
                    n += 1
                    if i + n == len(df_peak) - 1:
                        break
                    else:
                        pass
            elif df_peak[i] == 'STOP' and df_location[i] not in endlist:
                read_dens = 1
                moving_dens = 0.2
                n = 0
                move_operon = []
                move_operon.append(df_operon[i])
                while i > 0 and read_dens >= 0.2 and moving_dens >= 0.1:

                    index_start = full_location.index(df_location[i - n - 1])
                    index_end = full_location.index(df_location[i - n])
                    condition_values = values[index_start: index_end + 1]
                    if abs(df_location[i - n] - df_location[i - n - 1]) > 200:
                        for c in range(0, len(condition_values), 200):
                            for_condition_values = condition_values[c: c + 200]
                            non_zero_values = [z for z in for_condition_values if z != 0]
                            moving_dens = len(non_zero_values)/len(for_condition_values)
                            if moving_dens < 0.1:
                                break
                            else:
                                pass
                    else:
                        non_zero_values = [z for z in condition_values if z != 0]
                        read_dens = len(non_zero_values)/len(condition_values)


                    if read_dens >= 0.2 and moving_dens >= 0.1 and df_peak[i - n - 1] == 'START':
                        if '_a_half' in df_operon[i] or '_a_half' in df_operon[i - n - 1]:
                            if thestrand == 0:
                                mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df_location[i - n - 1]], 'end': [df_location[i]], 'score': [1000000], 'strand': ['+'], 'phase': ['.'], 'attributes': [name]}
                            else:
                                mydata = {'seqID': [genome_name], 'source': ['RendSeq'], 'feature': ['gene'], 'start': [df_location[i]], 'end': [df_location[i - n - 1]], 'score': [1000000], 'strand': ['-'], 'phase': ['.'], 'attributes': [name]}
                            newrow = pd.DataFrame(data = mydata)
                            GFF3_remain = pd.concat([GFF3_remain, newrow])
                            break
                        else:
                            move_operon.append(df_operon[i - n - 1])
                            if len(list(set(move_operon))) >= 2:
                                break
                            else:
                                pass
                    else:
                        pass
                    n += 1
                    if i - n == 0:
                        break
                    else:
                        pass
            else:
                pass

        GFF3 = pd.concat([GFF3, GFF3_remain])
        GFF3 = GFF3.drop_duplicates()
        GFF3 = GFF3.sort_values(['start'], axis = 0)
        GFF3 = GFF3.reset_index(drop = True)
        ### omit all annotations less than 15 nucleotides
        droplist = []
        startlist = GFF3['start'].tolist()
        endlist = GFF3['end'].tolist()
        i = 0
        for s, e in zip(startlist, endlist):
            if e - s < 15:
                droplist.append(i)
            else:
                pass
            i += 1
        GFF3 = GFF3.drop(GFF3.index[droplist])
        GFF3.to_csv(os.path.join(output_folder + os.path.basename(file) + '.gff3'), sep = '\t', index = False)

        del GFF3, GFF3_remain
    else:
        pass

if __name__ == '__main__':
    with Pool(8) as p:
        p.map(YourCode, sorted(glob.glob(input_folder + '*csv')))
