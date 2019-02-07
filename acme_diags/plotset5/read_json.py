import json

filename='/Users/zhang40/Documents/AIMS/repo/acme-diags/acme_diags/plotset5/obs_info_dictionary.json'
#filename='/Users/zhang40/Documents/AIMS/repo/acme-diags/acme_diags/plotset5/obs_info.json'

with open(filename) as json_data:
    data  = json.load(json_data)

print data['PRECT']['obs']
print data[0]
