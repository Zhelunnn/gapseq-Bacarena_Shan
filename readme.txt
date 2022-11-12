To get flux of metabolites:
1.The input is RDS file obtained from each cycle of Bacarena, Run scripts 'check_cf_0920.R' and 'crossfeeding_check_0919.py', which are programed by Shan. The fluxes of metabolite for each model can be checked in '*_sum_in_species.csv'.

To get fluxes of metabolites from reactions (extract reations fluxes from 'mflux'):
1.extract non-zero fluxes of reactions from 'mflux' and generate a two-column csv file (first column:reaction_id, second_column:flux)
2.taking csv file generated in first step as the input, run 'return_equation_bacarena_0709.py'. the output file contains all reactions with equation and flux. 