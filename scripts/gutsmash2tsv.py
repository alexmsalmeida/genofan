#!/usr/bin/env python

import os
import glob
import sys
from Bio import SeqIO

cluster_class_dict = {'Putrescine2spermidine': 'Aliphatic_amine', 'HGD_related': 'Putative', 'Arginine2putrescine': 'Aliphatic_amine', 'PFOR_II_pathway': 'SCFA', 'NADH_dehydrogenase_I': 'E-MGC', 'Respiratory_glycerol': 'E-MGC', 'Formate_dehydrogenase': 'E-MGC', 'Ech_complex': 'E-MGC', 'Nitrate_reductase': 'E-MGC','Molybdopterin_dependent_oxidoreductase':'E-MGC', 'Rnf_complex': 'E-MGC', 'Pyruvate2acetate-formate': 'SCFA','Sulfate2sulfide': 'E-MGC', 'Hydroxy-L-proline2proline': 'Other', 'Phenylacetate2toluene': 'Aromatic', 'Indoleacetate2scatole':'Aromatic', 'Fumarate2succinate': 'SCFA', 'Acetyl-CoA_pathway': 'SCFA','hydroxybenzoate2phenol': 'Aromatic', 'pdu': 'SCFA-other', 'EUT_pathway': 'SCFA-other', 'TMA': 'SCFA-other', 'p-cresol': 'Aromatic','Arginine2_Hcarbonate': 'Other', 'proline2aminovalerate': 'npAA', 'Leucine_reduction': 'SCFA', 'gallic_acid_met': 'Aromatic', 'bai_operon': 'Other', 'acetate2butyrate': 'SCFA', 'AAA_reductive_branch': 'Aromatic', 'porA': 'SCFA', 'Lysine_degradation': 'SCFA', 'glutamate2butyric': 'SCFA', 'caffeate_respiration': 'Aromatic', 'carnitine_degradaion_caiTABCDE': 'npAA', 'aminobutyrate2Butyrate': 'SCFA', 'succinate2propionate': 'SCFA', 'acrylate2propionate': 'SCFA', 'Threonine2propionate': 'SCFA', 'Glycine_reductase': 'SCFA', 'Glycine_cleavage': 'Other', 'histidine2glutamate_hutHGIU_operon': 'Other', 'succinate2propionate': 'SCFA', 'Oxidative_glycerol': 'Other', 'Flavoenzyme_AA_peptides_catabolism': 'Putative', 'Flavoenzyme_sugar_catabolism': 'Putative', 'Flavoenzyme_lipids_catabolism': 'Putative', 'OD_lactate_related': 'Putative', 'OD_eut_pdu_related': 'Putative', 'OD_AA_metabolism': 'Putative', 'OD_fatty_acids': 'Putative', 'OD_aldehydes_related': 'Putative', 'OD_unknown': 'Putative', 'TPP_fatty_acids': 'Putative', 'TPP_AA_metabolism': 'Putative', 'GR_AA_metabolism': 'Putative', 'GR_eut-pdu-related': 'Putative', 'GR_fatty_acids': 'Putative', 'OD_GR_eut_related': 'Putative', 'OD_GR_unassigned': 'Putative', 'Others_HGD_unassigned': 'Putative', 'fatty_acids-unassigned': 'Putative'}

genome_name = sys.argv[2]

for input_file in glob.glob(os.path.join(sys.argv[1], "*.gbk")):
    cluster_name = ".".join(os.path.basename(input_file).split(".")[:-1])
    with open(input_file) as file_in:
       for record in SeqIO.parse(file_in, "genbank"):
           for feature in record.features:
               if feature.type == "cand_cluster":
                   cluster_type = feature.qualifiers["product"]
                   cluster_class = []
                   for cluster in cluster_type:
                       cluster_class.append(cluster_class_dict[cluster])
                   print("%s\t%s\t%s\t%s" % (genome_name, cluster_name, "_".join(feature.qualifiers["product"]), "_".join(cluster_class)))
    
