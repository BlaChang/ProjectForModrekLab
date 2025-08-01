import json
import requests

main_categories = {
        "DNA Repair" : ["R-HSA-73894.5"],
        "DNA Replication" : ["R-HSA-69306.9"],
        "Apoptosis" : ["R-HSA-109581.6"],
        "Autophagy" : ["R-HSA-9612973.4"],
        "Cell Cycle Checkpoint": ["R-HSA-69620.5"],
        "Cellular Response to Stress": ["R-HSA-2262752.13"],
        "Cell Cell Communication" : ["R-HSA-1500931.8"],
        "Cytokine Signalling Immune System" : ["R-HSA-1280215.7"],
        "Gene Expression(Transcription)": ["R-HSA-74160.9"],
        "Hypoxia" : ["R-HSA-1234174.5"],
        "Metabolism of Carbs" : ["R-HSA-71387.14"],
        "Mitochondrial Stress Response" : ["R-HSA-9840373.2"],
        "Metabolism of lipids" : ["R-HSA-556833.9"],
        "Metabolism of Proteins" : ["R-HSA-392499.12"],
        "Metabolism of Amino Acids and Derivatives" : ["R-HSA-71291.10"],
        "Metabolism of RNA" : ["R-HSA-8953854.10"],
        "Epigenetics" : ["R-HSA-4839726.5", "R-HSA-212165.7"],
        "Stemness" : ["R-HSA-452838.7"],
        "Proliferation" : ["R-HSA-2892247.5"],
        "Differentiation" : ["R-HSA-9758941.6"],
        "Extracellular Matrix Formation" : ["R-HSA-1474244.5"],
        "Adaptive Immune System" : ["R-HSA-1280218.8"],
        "Innate Immune System" : ["R-HSA-168249.10"],
        "Angiogenesis" : ["R-HSA-210993.3"],
        "Protein Localization" : ["R-HSA-9609507.4"],
        "Signal Transduction" : ["R-HSA-162582.13"]
        }
category_dict = dict()
for category in main_categories.keys():
    pathways = []
    for reactome_id in main_categories[category]:
        request_string = f'https://reactome.org/ContentService/data/pathway/{reactome_id}/containedEvents/displayName'
        headers = {
                'User-Agent': 'python-requests/2.31.0',
                'Accept-Encoding': 'gzip, deflate',
                'Accept': '*/*',
                'Connection': 'keep-alive',
                }
        response = requests.get(request_string, headers=headers)
        contained_pathways = response.text[1:-1].split(",")
        contained_pathways = [x.strip() for x in contained_pathways]
        pathways = pathways + contained_pathways
    category_dict[category] = pathways


with open("more_categories.json", "w") as f:
    json.dump(category_dict,f, indent=4)
    
