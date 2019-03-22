#! /usr/bin/env bash

## format the kegg file so we can use it
python3 kegg_database_format.py hsa00001.keg

## get the hit lists
# hit_maker.py results_list.txt name_or_ID pvalue
python3 hit_maker.py /Users/Jacob/Box/recurrent\ recurrent\ positive\ selection/190312_results/2019-03-05-10-40_Hsap_results.parsed.txt name 3.971721344030503e-06
python3 hit_maker.py /Users/Jacob/Box/recurrent\ recurrent\ positive\ selection/190312_results/2019-03-05-10-55_Hsap9_results.parsed.txt name 4.064379775646237e-06
python3 hit_maker.py /Users/Jacob/Box/recurrent\ recurrent\ positive\ selection/190312_results/2019-03-05-11-08_Mmus_results.parsed.txt ID 3.17017499365965e-06
python3 hit_maker.py /Users/Jacob/Box/recurrent\ recurrent\ positive\ selection/190312_results/2019-03-05-12-12_Pman_results.parsed.txt ID 3.631345776744862e-06
python3 hit_maker.py /Users/Jacob/Box/recurrent\ recurrent\ positive\ selection/190312_results/2019-03-05-12-28_Mluc_results.parsed.txt ID 4.305148958153952e-06
python3 hit_maker.py /Users/Jacob/Box/recurrent\ recurrent\ positive\ selection/190312_results/2019-03-05-12-49_Btau_results.parsed.txt ID 1.5422578655151142e-05
python3 hit_maker.py /Users/Jacob/Box/recurrent\ recurrent\ positive\ selection/190312_results/2019-03-05-12-51_Clup_results.parsed.txt ID 7.627765064836003e-06

## convert all the hit lists into human gene names
# ortho_conversion.py hit_list.txt orthology_table.txt human_ensembles.csv
python3 ortho_conversion.py /Users/Jacob/corsair/ontology_analysis/hit_lists/Btau_hits.txt /Users/Jacob/corsair/ontology_analysis/conversion_tables/human_to_cow.txt /Users/Jacob/corsair/ontology_analysis/human_cds_ensemble.csv
python3 ortho_conversion.py /Users/Jacob/corsair/ontology_analysis/hit_lists/Clup_hits.txt /Users/Jacob/corsair/ontology_analysis/conversion_tables/human_to_dog.txt /Users/Jacob/corsair/ontology_analysis/human_cds_ensemble.csv
python3 ortho_conversion.py /Users/Jacob/corsair/ontology_analysis/hit_lists/Mluc_hits.txt /Users/Jacob/corsair/ontology_analysis/conversion_tables/human_to_microbat.txt /Users/Jacob/corsair/ontology_analysis/human_cds_ensemble.csv
python3 ortho_conversion.py /Users/Jacob/corsair/ontology_analysis/hit_lists/Mmus_hits.txt /Users/Jacob/corsair/ontology_analysis/conversion_tables/human_to_mouse.txt /Users/Jacob/corsair/ontology_analysis/human_cds_ensemble.csv
python3 ortho_conversion.py /Users/Jacob/corsair/ontology_analysis/hit_lists/Pman_hits.txt /Users/Jacob/corsair/ontology_analysis/conversion_tables/human_to_deermouse.txt /Users/Jacob/corsair/ontology_analysis/human_cds_ensemble.csv

## get the kegg pathway counts for each clade
# multi_dics.py converted_hit_list.txt
python3 multi_dics.py /Users/Jacob/corsair/ontology_analysis/converted_hit_lists/Hsap_converted_hits.txt
python3 multi_dics.py /Users/Jacob/corsair/ontology_analysis/converted_hit_lists/Mluc_converted_hits.txt
python3 multi_dics.py /Users/Jacob/corsair/ontology_analysis/converted_hit_lists/Pman_converted_hits.txt
python3 multi_dics.py /Users/Jacob/corsair/ontology_analysis/converted_hit_lists/Mmus_converted_hits.txt
python3 multi_dics.py /Users/Jacob/corsair/ontology_analysis/converted_hit_lists/Btau_converted_hits.txt
python3 multi_dics.py /Users/Jacob/corsair/ontology_analysis/converted_hit_lists/Clup_converted_hits.txt
