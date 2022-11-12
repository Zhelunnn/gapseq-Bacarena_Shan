import argparse
#################add equation information for the rxn file extracted from bacarena model#############
##########rxn file was the reactions that have the unzero flux####
def read_equation_from_db_rxn(eqa):#######read fifth column of seed_reactions_corrected file
    cpd_list=eqa.split(';')
    consume_list=[] # in form: [-1cpd00001[c0],-2cpd00002[c0]]
    product_list=[]    # [ 1cpd00003[e0],3cpd00004[c0]]
    for cpd in cpd_list:
        cpd_coeff=cpd.split(':')[0]
        cpd_name=cpd.split(':')[1]
        if cpd.split(':')[2]=='0':
            cpd_loc='[c0]'
        elif cpd.split(':')[2]=='1':
            cpd_loc='[e0]'
        else:
            cpd_loc='[p0]'
        coe_cpd=''.join([cpd_coeff,cpd_name,cpd_loc])
        if cpd.startswith('-'):
            consume_list.append(coe_cpd[1:])
        else:
            product_list.append(coe_cpd)
    return consume_list,product_list

###########generate active file from bacarena model###########
#example:rxn00173_c0\tacetyl-CoA:phosphate acetyltransferase\t1cpd00009[c0]+1cpd00022[c0]>1cpd00010[c0]+1cpd00196[c0]\t1000
def equation_re_arena(bac_file,db_file,output_file):
    new_file=open(output_file,'w')
    bac_rxn_list=open(bac_file).readlines()#    bac file in form: rxn_id\tflux\n
    for rxn_line in bac_rxn_list[1:]:#exclude the first blank line
        rxn_flux = rxn_line.strip('\n').split(',')[1]
        if float(rxn_flux)!=0:
            rxn_id=rxn_line.strip('\n').split(',')[0][1:9]
            for db_rxn_line in open(db_file):
                if rxn_id in db_rxn_line.split('\t')[0]:##avoid the appearance of rxnid in other columns
                    eqa_original=db_rxn_line.split('\t')[4]
                    eqa_name=db_rxn_line.split('\t')[2]
                    consu='+'.join(read_equation_from_db_rxn(eqa_original)[0])
                    prod='+'.join(read_equation_from_db_rxn(eqa_original)[1])
                    new_file.write('%s\t%s\t%s>%s\t%s\n' %(rxn_id,eqa_name,consu,prod,rxn_flux))
                    break
    new_file.close()
# if __name__=='__main__':
#     active_reactions=argparse.ArgumentParser()
#     active_reactions.add_argument('-i',required=True,help='input file should be grenerated from bacarena model')
#     active_reactions.add_argument('-o',required=False,default='bac_rxns.txt',help='output file was required and the format is rxn,coeffi+cpdid,fluxes')
#     active_reactions.add_argument('-db',required=False,default='/srv/scratch/z5245780/software/gapseq/gapseq_1.2/dat/seed_reactions_corrected.tsv',help='the database used to add the rxn equation')
#     args = vars(active_reactions.parse_args())
#     equation_re_arena(args['i'],args['db'],args['o'])

#equation_re_arena('reaction_2.tsv','seed_reactions_corrected.tsv','123.txt')