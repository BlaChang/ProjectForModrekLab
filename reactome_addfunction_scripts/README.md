#Add main function column(and common name column) to tsv
python add_main_function.py INPUTFILE

#generate categories json
python generate_categoies.json

#generate msig-reactome.json
python reactome_name_dict.py
