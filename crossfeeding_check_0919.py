import argparse
import os
parser = argparse.ArgumentParser()
parser.add_argument('-E.g',
                    required=False,
                    help='python3 /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Get_network_graph/Get_network_graph_Willis.py -indir -infile')
parser.add_argument('-mydir',default='/srv/scratch/z5245780/SpongeMAGs_Robbins_ISME/all_sponge_DB_bins_over50pctComp/CAR/genome/Bacarena_community',
                    required=False,
                    help='Full path of input files.')
parser.add_argument('-infile',
                    required=True,
                    help='Name of input files.(generated from RDS by the R script')

args = vars(parser.parse_args())

mydir = args['mydir']
infile_original = args['infile']

#################### Or run script within pycharm:
######################## user define starts #########################
# mydir = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/Github_repsitory/Shan-George-master/BacArena/BacArena_STY_v20210520_3autotrophs_OTU2n6n8_400grids_mineral_sw/Katana/R_3.6.3__BacArena_1.8.2_CRAN/mineral_sw_210716_8taxa_BacArena_simplify_ExImf_false_fba/Rscript_BacArena__Best_models_STY_8taxa_20210902__mtf_rmRate_tp_ExInf_replenish_growlimit_parL_diffModels_addSubsF_auto_2/8taxa__v210903_sw_HT_Be_915_3__NH3500_HT30_rplc2/'
# infile_original = '100_exchangelist_HT_Be_915_3__NH3500_HT30_rplc2_rl_50_cpds60_time2.csv'
######################## user define ends #########################
outfile_1_cpds_sum_in_species = mydir + '/' + infile_original.split('.csv')[0] + '_sum_in_species.csv'
outfile_2_network_format      = mydir + '/' + infile_original.split('.csv')[0] + '_network.csv'

#
# st1. get 8 dicts of dicts: for the first dict, each species are keys,
# for the second dict, the cpd IDs are the keys and fluxes are values.
# At last, by providing species, will get you  its flux of each metabolites.
dict_of_dict ={}
header_list = []
#for each_line in open(mydir +'/'+ infile_original):
for each_line in open(mydir + '/'+infile_original):
    each_line_split = each_line.strip('\n').split(',')

    if each_line.startswith('"species",') or each_line.startswith('species,') :
        header_list = each_line_split
    else:
        species_id = each_line_split[0]

        if species_id not in dict_of_dict:
            dict_of_dict[species_id] = {}

        for (cpd, flx) in zip(header_list[1:], each_line_split[1:]):
            if flx != 'NA':

                if cpd not in dict_of_dict[species_id]:
                    dict_of_dict[species_id][cpd] = float(flx)
                else:
                    dict_of_dict[species_id][cpd] += float(flx)
                    
# st2. write the table of sum of different cpds in each species to an output file.
outfile_handle = open(outfile_1_cpds_sum_in_species, 'w')
outfile_handle.write(','.join(header_list) + '\n')
for each_species in dict_of_dict:
    flx_sum_list = []
    for each_cpd in header_list[1:]:
        flx_sum_list.append(dict_of_dict[each_species].get(each_cpd, 'NA'))

    flx_sum_list_str = [str(i) for i in flx_sum_list]


    outfile_handle.write('%s,%s\n' % (each_species, ','.join(flx_sum_list_str)))
outfile_handle.close()


# st3. write a table in the form for network graph in r.
outfile_handle = open(outfile_2_network_format, 'w')
outfile_handle.write('%s,%s,%s,%s,%s\n' % ('met', 'prod', 'cons', 'prod.flux', 'con.flux'))

for each_cpd in header_list:

    provide_list = []
    comsumer_list = []
    for each_species in dict_of_dict:
        # print('%s\t%s' % (each_species, dict_of_dict[each_species].get(each_cpd, 0)))
        flx_sum = dict_of_dict[each_species].get(each_cpd, 0)
        if flx_sum > 0:
            provide_list.append(each_species)
        if flx_sum < 0:
            comsumer_list.append(each_species)

    if (provide_list == []) and (comsumer_list != []):
        provide_list = ["Environment"]
        # print('%s\t%s\t%s' % (each_cpd, provide_list, comsumer_list))
    if (provide_list != []) and (comsumer_list == []):
        comsumer_list = ["Environment"]
        # print('%s\t%s\t%s' % (each_cpd, provide_list, comsumer_list))

    if (provide_list != []) and (comsumer_list != []):
        # pass
        #print('%s\t%s\t%s' % (each_cpd, provide_list, comsumer_list))

        for each_pro in provide_list:
            for each_con in comsumer_list:

                pro_flux = dict_of_dict.get(each_pro, {}).get(each_cpd, 0) # provide blanks if keys can not be found in the dictionaries.
                con_flux = dict_of_dict.get(each_con, {}).get(each_cpd, 0)

                outfile_handle.write('%s,%s,%s,%s,%s\n' % (each_cpd, each_pro, each_con, pro_flux, con_flux))

outfile_handle.close()

###########code from here to the top can be used to check the specific cycle crossfeeding, need to add the species at the beginning of first line and modify 'prod_flux' to 'prod.flux'...#######
##########change cpd_id to met, providers to prod, consumers to cons.#########
##################### Or run script within pycharm:
######################### user define starts #########################

dict_seed_metabolites_edited_gapseq = '/srv/scratch/z5245780/software/gapseq/gapseq_1.2/dat/seed_metabolites_edited.tsv'

######################### user define ends #########################
infile_2  = mydir + '/' + infile_original.split('.csv')[0] + '_network.csv'
outfile_3 = mydir + '/' + infile_original.split('.csv')[0] + '_out_final_cf.csv'

# st1. create a dict of MS metabolite id to metabolite name.
dict_metid_2_metnm = {}
for each in open(dict_seed_metabolites_edited_gapseq):
    each_split = each.strip().split('\t')
    if not each.startswith('id'):
        met_id = each_split[0]
        met_nm = each_split[3]
        dict_metid_2_metnm[met_id] = met_nm
        # print(met_id)

# st2. use met_id to get met_nm from the input file.
outfile_handle = open(outfile_3, 'w')
for each in open(infile_2):
    each_split = each.strip().split(',')
    if each.startswith('met'):
        outfile_handle.write('%s,%s,%s\n' % (each_split[0], 'cpd_nm', ','.join(each_split[1:])))
    else:
        met_id = each_split[0].split('EX_')[1].split('_')[0]
        # print(met_id)
        met_nm = dict_metid_2_metnm[met_id]
        outfile_handle.write('%s,%s,%s\n' % (met_id, met_nm, ','.join(each_split[1:])))
outfile_handle.close()
os.system('rm '+infile_2)
# print('Then run the script /Users/zzfanyi/PycharmProjects/EMP500/GapSeq/Replace_modelID_with_Species_nm/Replace_string__modelID_with_Species_nm.py to replace model IDs with species names.')

#
# # @2021-09-24
# ##################### Or run script within pycharm:
# ######################### user define starts #########################
# # dict = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/Dict_STY_8taxa_modelID_species_nm.txt'
# # mydir = '/Users/zzfanyi/Bioinfo/EMP500/02-shotgun/GapSeq_models/Github_repsitory/Shan-George-master/BacArena/BacArena_STY_v20210520_3autotrophs_OTU2n6n8_400grids_mineral_sw/Katana/R_3.6.3__BacArena_1.8.2_CRAN/mineral_sw_210716_8taxa_BacArena_simplify_ExImf_false_fba/Rscript_BacArena__Best_models_STY_8taxa_20210902__mtf_rmRate_tp_ExInf_replenish_growlimit_parL_diffModels_addSubsF_auto_2/8taxa__v210903_sw_HT_Be_915_3__NH3500_HT30_rplc2/'
# # infile_original = '100_exchangelist_HT_Be_915_3__NH3500_HT30_rplc2_rl_50_cpds60_time2.csv'
# ######################### user define ends #########################
# infile_3  = mydir + '/' + infile_original.split('.csv')[0] + '_out3.csv'
# outfile_4 = mydir + '/' + infile_original.split('.csv')[0] + '_out4_network_form.csv'
#
# outfile_handle = open(outfile_4, 'w')
# for each in open(infile_3):
#     # print(each)
#     each_new1 = each.replace('"STY_Merged_OTU01"', 'Ca. Spongiicolonus versatilis')
#     each_new2 = each_new1.replace('"STY_Merged_OTU02"', 'Ca. Synechococcus stylissus')
#     each_new3 = each_new2.replace('"STY_Merged_OTU03"', 'Ca. Deltaporibacter stylissus')
#     each_new4 = each_new3.replace('"STY_Merged_OTU04"', 'Ca. Gammaporibacter thiotrophicus')
#     each_new5 = each_new4.replace('"STY_Merged_OTU05"', 'Ca. Oxydemutator thiolithooxidans')
#     each_new6 = each_new5.replace('"STY_Merged_OTU06"', 'Ca. Nitrospongiibacter stylissus')
#     each_new7 = each_new6.replace('"STY_Merged_OTU07"', 'Ca. Spongiihabitans thiooxidans')
#     each_new8 = each_new7.replace('"STY_Merged_OTU08"', 'Ca. Cenoporarchaeum stylissus')
#
#     each_new8 = each_new8.replace('STY_Merged_OTU01', 'Ca. Spongiicolonus versatilis')
#     each_new8 = each_new8.replace('STY_Merged_OTU02', 'Ca. Synechococcus stylissus')
#     each_new8 = each_new8.replace('STY_Merged_OTU03', 'Ca. Deltaporibacter stylissus')
#     each_new8 = each_new8.replace('STY_Merged_OTU04', 'Ca. Gammaporibacter thiotrophicus')
#     each_new8 = each_new8.replace('STY_Merged_OTU05', 'Ca. Oxydemutator thiolithooxidans')
#     each_new8 = each_new8.replace('STY_Merged_OTU06', 'Ca. Nitrospongiibacter stylissus')
#     each_new8 = each_new8.replace('STY_Merged_OTU07', 'Ca. Spongiihabitans thiooxidans')
#     each_new8 = each_new8.replace('STY_Merged_OTU08', 'Ca. Cenoporarchaeum stylissus')
#
#     # print(each_new8)
#     outfile_handle.write(each_new8)
# outfile_handle.close()
# print("Done! Now enjoy make your network graphs! :D")
# print("Use ***_out4_network_form.csv for your next step~")