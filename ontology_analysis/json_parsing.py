#!/usr/bin/env python3

import json
import sys
import copy

json_file = sys.argv[1]

with open(json_file, 'r') as f:
    kegg = json.load(f)

# kegg_new = {}
# for a in kegg['children']:
#     kegg_new[a['name']] = a['children']
# kegg = copy.deepcopy(kegg_new)
# kegg_new = {}

# for a in kegg:
#     kegg_new[a] = {}
#     for dic in kegg[a]:
#         kegg_new[a][dic['name']] = dic['children']
# kegg = copy.deepcopy(kegg_new)
# kegg_new = {}

# for a in kegg:
#     kegg_new[a] = {}
#     for b in kegg[a]:
#         kegg_new[a][b] = {}
#         for dic in kegg[a][b]:
#             if 'children' in dic:
#                 kegg_new[a][b][dic['name']] = dic['children']
#             else:
#                 kegg_new[a][b][dic['name']] = {}
# print(kegg_new)




def reduce_json(dic):
    new_dict = {}
    
    if ("name" in dic) and ("children" in dic):
        new_dict[dic['name']] = {}
        for next_level in dic['children']:
            new_dict[dic['name']][next_level['name']] = reduce_json(dic['children'])
    else:
        return(dic['name'])
    return new_dict

# print(kegg.keys())

kegg = reduce_json(kegg)
for key in kegg:
    print(kegg[key].keys())