import pandas as pd
import json
import subprocess
import re
import argparse

def add_commonname_column(input_path, output_path):
    df = pd.read_csv(input_path, sep = "\t")
    with open("msig-reactome.json","r") as f:
        with open("more_categories.json","r") as g:
            reactome_dict = json.load(f)   
            categories = json.load(g)

    #Custom Columns need to be added due to being in an archived/outdated gene setj
    reactome_dict["REACTOME_CD28_CO_STIMULATION"] = "Co-stimulation by CD28"
    reactome_dict["REACTOME_CHROMATIN_MODIFYING_ENZYMES"] = "Chromatin modifying enzymes"
    reactome_dict["REACTOME_COSTIMULATION_BY_THE_CD28_FAMILY"] = "Regulation of T cell activation by CD28 family"
    reactome_dict["REACTOME_DEVELOPMENTAL_BIOLOGY"] = "Developmental Biology"
    reactome_dict["REACTOME_MATURATION_OF_SPIKE_PROTEIN"] = "Maturation of Spike Protein"
    reactome_dict["REACTOME_RORA_ACTIVATES_GENE_EXPRESSION"] = "None"
    reactome_dict["REACTOME_CTLA4_INHIBITORY_SIGNALING"] = "Co-inhibition by CTLA4"
    reactome_dict["REACTOME_PD_1_SIGNALING"] = "Co-inhibition by PD-1"
    reactome_dict["REACTOME_CARNITINE_METABOLISM"] = "Carnitine shuttle"
    reactome_dict["REACTOME_METABOLISM_OF_CARBOHYDRATES"] = " Metabolism of carbohydrates and carbohydrate derivatives"
    reactome_dict["REACTOME_SERINE_BIOSYNTHESIS"] = "Serine metabolism"
    reactome_dict["REACTOME_BMAL1_CLOCK_NPAS2_ACTIVATES_CIRCADIAN_GENE_EXPRESSION"] = " BMAL1:CLOCK,NPAS2 activates circadian expression"
    reactome_dict["REACTOME_ARACHIDONIC_ACID_METABOLISM"] = "Arachidonate metabolism"
    
    #More Custom Columns that I will add except this is covers the outdated genesets for EVERY pathway, common, unique, reallyunique, semiunique whatever lol..
    reactome_dict["REACTOME_ACTIVATION_OF_IRF3_IRF7_MEDIATED_BY_TBK1_IKK_EPSILON"] = "Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE)"
    reactome_dict["REACTOME_CDT1_ASSOCIATION_WITH_THE_CDC6_ORC_ORIGIN_COMPLEX"] = 'CDT1 association with the CDC6:ORC:origin complex'
    reactome_dict["REACTOME_DEATH_RECEPTOR_SIGNALLING"] = "Death Receptor Signaling"
    reactome_dict["REACTOME_DEPOLYMERISATION_OF_THE_NUCLEAR_LAMINA"] = "Depolymerization of the Nuclear Lamina"
    reactome_dict["REACTOME_DSCAM_INTERACTIONS"] = 'DSCAM interactions'
    reactome_dict["REACTOME_GAMMA_CARBOXYLATION_HYPUSINE_FORMATION_AND_ARYLSULFATASE_ACTIVATION"] = "Gamma carboxylation, hypusinylation, hydroxylation, and arylsulfatase activation"
    reactome_dict["REACTOME_GLYCOLYSIS"] = 'Glycolysis'
    reactome_dict["REACTOME_NICOTINAMIDE_SALVAGING"] = "Don't exist anymore"
    reactome_dict["REACTOME_PYRUVATE_METABOLISM_AND_CITRIC_ACID_TCA_CYCLE"] = 'Pyruvate metabolism and Citric Acid (TCA) cycle'
    reactome_dict["REACTOME_REGULATION_OF_TNFR1_SIGNALING"] = "Regulation of TNFR1 signaling"
    reactome_dict["REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS"] = 'Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.'
    reactome_dict["REACTOME_SIGNALING_BY_FGFR3_FUSIONS_IN_CANCER"] = "Signaling by FGFR3 fusions in cancer"
    reactome_dict["REACTOME_TAK1_ACTIVATES_NFKB_BY_PHOSPHORYLATION_AND_ACTIVATION_OF_IKKS_COMPLEX"] = "TAK1 activates NFkB by phosphorylation and activation of IKKs complex"
    reactome_dict["REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT"] = "Aerobic respiration and respiratory electron transport"
    reactome_dict["REACTOME_TNFR1_INDUCED_NFKAPPAB_SIGNALING_PATHWAY"] = "TNFR1-induced NF-kappa-B signaling pathway"
    main_function = []
    main_function_two = []
    commonname = []
    for pathway in df["pathway"]:
        i = 0
        for category in categories.keys():
            try:
                common_name = reactome_dict[pathway]
            except:
                try:
                    #Hacky script to try and get more reactome names, when the gene set is outdated, this sometimes works.
                    url = f"https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/{pathway}.html" 
                    reactome_id_unfiltered = subprocess.run(f"curl  {url} | grep Brief -A1", capture_output = True, shell = True, text = True)
                    regex_capture = r"(?s:.)+\<td\>(.+)\<\/td(?s:.)+"
                    match = re.search(regex_capture, reactome_id_unfiltered.stdout)
                    reactome_dict[pathway] = match.group(1)
                    print(f"""reactome_dict["{pathway}"] = '{match.group(1)}'""")
                except:
                    print(f"""reactome_dict["{pathway}"] = """)
                    reactome_dict[pathway] = "Doesn't work"
                    break
            if common_name in categories[category]:
                i+=1
                if i == 1:
                    main_function.append(category)
                    commonname.append(reactome_dict[pathway])
                elif i == 2:
                    print("Multiple main_function detected")
                    main_function_two.append(category)
        if i ==0:
            main_function.append("Other")
            main_function_two.append("None")
            try:
                commonname.append(reactome_dict[pathway])
            except:
                commonname.append("failed")
        elif i==1:
            main_function_two.append("None")

    df.insert(1, "main_function", main_function)
    df.insert(0, "common_name", commonname)
    df.insert(3,"main_function_two",main_function_two)
    df.drop(axis = 0, index = 0, inplace = True)
    df.to_csv(output_path, sep = "\t", index_label = False, index = False)
    print(f"Saved to {output_path}")
    print("###FIRST FIVE COLUMNS")
    print(df.head())



def main():
    parser = argparse.ArgumentParser(description="Just an example",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input")
    parser.add_argument("-o", "--output", default = "output.txt")
    args = parser.parse_args()
    args_dict = vars(args)
    if args_dict["output"] == "output.txt":
        args_dict["output"] = args_dict["input"].split(".")[0] + "_with_mainfunction.txt"
    add_commonname_column(args_dict["input"], args_dict["output"])

    

if __name__ == '__main__':
    main()




