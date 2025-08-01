import json
import requests

# Open and read the JSON file
with open('c2.cp.reactome.v2025.1.Hs.json', 'r') as file:
    data = json.load(file)

pretty_json_string = json.dumps(data, indent=4)
names = data.keys()
reactome_dict = dict()
for name in names:
    reactome_id = data[name]["exactSource"]
    url = f"https://reactome.org/ContentService/data/query/{reactome_id}/displayName"
    headers = {
            'User-Agent': 'python-requests/2.31.0',
            'Accept-Encoding': 'gzip, deflate',
            'Accept': '*/*',
            'Connection': 'keep-alive',
            }
    response = requests.get(url, headers=headers)
    reactome_dict[name] = response.text
    reactome_dict[response.text] = name
    print(response.text)
    with open("data.json", "w") as f:
        json.dump(reactome_dict,f, indent=4)

print(reactome_dict)
