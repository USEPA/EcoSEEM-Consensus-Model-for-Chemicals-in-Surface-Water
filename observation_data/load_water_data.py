# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 22:16:34 2020
#https://netbeans.org/kb/docs/ide/install-and-configure-mysql-server.html

A. Load data from National Water Quality Data Portal
A.1. Load and filter station data


@author: Risa Sayre
"""
import numpy as np
#import time
import pandas as pd
#import mysql.connector
#from sqlalchemy import create_engine
#import sys
import os
import xml.etree.ElementTree as ET
import pickle as pkl

os.chdir('C:\\Users\\rrace\\Documents\\ECOSEEM')

# =============================================================================
# 
# # do once: find unique chem names
# start_year = 2008
# end_year = 2019
# date_range = [i for i in range(start_year,end_year)]
# 
# df = pd.read_csv('narrowresult_2008.zip', compression='zip', sep=',')
# val_counts = df.CharacteristicName.value_counts()
# chem_count = (pd.DataFrame({'CharacteristicName': val_counts.index,
#                             '2008': val_counts.values}))
# 
# for d in date_range:
#     filename = 'narrowresult_'+ str(d) +'.zip'
#     df = pd.read_csv(filename, compression='zip', sep=',', low_memory=False)
#     val_counts = df.CharacteristicName.value_counts()
#     char_name_val_count = (pd.DataFrame({'CharacteristicName': val_counts.index,
#                                          d: val_counts.values}))
#     chem_count = pd.merge(chem_count,
#                           char_name_val_count,
#                           on='CharacteristicName',
#                           how='outer')
# =============================================================================

# =============================================================================
# # do once: map to dtxsid and sum to find high frequency names
# # from Dashboard, batch download TEST, OPERA, QSAR-ready SMILES, average mass
# dtxsids = pd.read_excel('dtxsids_from_dashboard.xls')
# result = pd.merge(chem_count,
#                   dtxsids,
#                   how='left',
#                   left_on='CharacteristicName',
#                   right_on='INPUT')
# =============================================================================


bad_CharacteristicNames = ['Oil and Grease','Oil and grease',
                           'Methylmercury(1+)','C1-C4 Fluoranthenes',
                           'Trihalomethanes','Petroleum phenols', 
                           'Phenols','Total Sulfate',
                           'C1-C4 Phenanthrenes',
                           'Polychlorinated biphenyls','C1-C4 Chrysenes',
                           'Benzene, 1,2,4,5-tetrachloro- and/or 1,2,3,5-Tetrachlorobenzene',
                           'Cylindrospermopsins','Hydrocarbons, petroleum',
                           'Aroclor 1016 mixt. with Aroclor 1242',
                           'Total microcystins plus nodularins',
                           'Diesel range organics','Chlordane, technical, and/or chlordane metabolites',
                           'Tannin and Lignin','Organics volatile mix, unspecified', 
                           'Hexazinone Transformation Product C','Gasoline range organics','Diuron Metabolite',
                           'Hexazinone Transformation Product D',
                           "DDT/DDD/DDE, sum of p,p' & o,p' isomers",
                           "DDT/DDD/DDE, sum of p,p' isomers",
                           'Chlorthal monoacid and diacid degradates',
                           'Transphilic fraction of organic carbon','Specific UV Absorbance at 254 nm',
                           'BZ#001L',
                           'BZ#003L',
                           'BZ#004L',
                           'BZ#015L',
                           'BZ#019L',
                           'BZ#028L',
                           'BZ#037L',
                           'BZ#054L',
                           'BZ#077L',
                           'BZ#081L',
                           'BZ#104L',
                           'BZ#105L',
                           'BZ#111L',
                           'BZ#114L',
                           'BZ#118L',
                           'BZ#123L',
                           'BZ#126L',
                           'BZ#155L',
                           'BZ#156L',
                           'BZ#157L',
                           'BZ#167L',
                           'BZ#169L',
                           'BZ#178L',
                           'BZ#188L',
                           'BZ#189L',
                           'BZ#202L',
                           'BZ#205L',
                           'BZ#206L',
                           'BZ#208L',
                           'BZ#209L',                           
                           'Hexazinone Transformation Product E','Sulfur Sulfate',
                           #mixtures
                           'ortho & para Xylene mix',
                           "Surfactants, unspecified mix","Phosphated pesticides", "Phenols",
                           "Xylenes, m- & p- Mix***retired***use m,p-Xylene",
                           'meta & para Xylene mix***retired***use m,p-Xylene',"Saxitoxins",'Aliphatics fraction',
                           'Aroclor (unspecified)',
                           'Aroclor 1016 mixt. with Aroclor 1221',
                           'Aroclor 1016 mixt. with Aroclor 1242',
                           'Aromatics fraction',
                           'Aluminum, Organic + Inorganic Monomeric (reactive aluminum)',
                           'Aluminum, Organic Monomeric (reactive aluminum)',
                           'Colored dissolved organic matter (CDOM)',
                           'Cyhalothrins','4-Nonylphenols',
                           'Hydrophilic fraction of organic carbon',
                           'Haloacetic acids',
                           'Halogenated organics',
                           'Gasoline range organics',
                           'Dissolved microcystins plus nodularins',
                           'Diesel range organics',
                           'Hydrophobic fraction of organic carbon',
                           'Hexazinone Transformation Product C',
                           'Hexazinone Transformation Product D',
                           'Hexazinone Transformation Product E',
                           'Hydrocarbons, C>28-C40',
                           'Hydrocarbons, C10-28',
                           'Hydrocarbons, petroleum',
                           'Hydrocarbons, total purgeable',
                           'Hydrocarbons, volatile',
                           'o-Xylene, mixt. with m-xylene and p-xylene',
                           'PCB-105/132/153',
                           'PCB-106/118',
                           'PCB-107/108',
                           'PCB-107/109',
                           'PCB-107/124',
                           'PCB-108/112',
                           'PCB-110/115',
                           'PCB-111/115',
                           'PCB-112/119',
                           'PCB-114/131',
                           'PCB-115/87',
                           'PCB-12/13',
                           'PCB-122/131',
                           'PCB-123/149',
                           'PCB-124/135/144',
                           'PCB-124/147',
                           'PCB-128/162',
                           'PCB-128/166',
                           'PCB-128/167',
                           'PCB-129/138/160/163',
                           'PCB-129/138/163',
                           'PCB-129/178',
                           'PCB-131/165',
                           'PCB-132/153',
                           'PCB-132/161',
                           'PCB-132/168',
                           'PCB-133/142',
                           'PCB-134/143',
                           'PCB-135/144',
                           'PCB-135/151',
                           'PCB-135/151/154',
                           'PCB-137/176',
                           'PCB-138/158',
                           'PCB-138/163',
                           'PCB-138/163/164',
                           'PCB-139/140',
                           'PCB-139/149',
                           'PCB-146/165',
                           'PCB-147/149',
                           'PCB-15/17',
                           'PCB-153/168',
                           'PCB-156/157',
                           'PCB-156/171/202',
                           'PCB-157/200',
                           'PCB-158/160',
                           'PCB-16/32',
                           'PCB-163/164',
                           'PCB-170/190',
                           'PCB-171/173',
                           'PCB-171/202',
                           'PCB-172/197',
                           'PCB-18/30',
                           'PCB-180/193',
                           'PCB-182/187',
                           'PCB-183/185',
                           'PCB-194/205',
                           'PCB-195/208',
                           'PCB-196/203',
                           'PCB-197/200',
                           'PCB-198/199',
                           'PCB-20/21/33',
                           'PCB-20/28',
                           'PCB-20/33',
                           'PCB-21/33',
                           'PCB-21/33/53',
                           'PCB-24/27',
                           'PCB-26/29',
                           'PCB-27/33/53',
                           'PCB-28/31',
                           'PCB-37/42',
                           'PCB-4/10',
                           'PCB-40/41/71',
                           'PCB-41/64',
                           'PCB-41/64/71',
                           'PCB-41/64/71/72',
                           'PCB-41/71',
                           'PCB-42/59',
                           'PCB-43/49',
                           'PCB-43/52',
                           'PCB-44/47/65',
                           'PCB-45/51',
                           'PCB-47/48',
                           'PCB-47/48/75',
                           'PCB-48/75',
                           'PCB-49/69',
                           'PCB-5/8',
                           'PCB-50/53',
                           'PCB-52/69',
                           'PCB-56/60',
                           'PCB-59/62/75',
                           'PCB-61/70',
                           'PCB-61/70/74/76',
                           'PCB-61/74',
                           'PCB-65/75',
                           'PCB-66/76',
                           'PCB-66/95',
                           'PCB-7/9',
                           'PCB-70/76',
                           'PCB-77/110',
                           'PCB-81/87',
                           'PCB-82/151',
                           'PCB-83/99',
                           'PCB-84/92',
                           'PCB-85/116',
                           'PCB-85/116/117',
                           'PCB-86/87/97/108/119/125',
                           'PCB-87/115',
                           'PCB-87/117/125',
                           'PCB-88/91',
                           'PCB-90/101',
                           'PCB-90/101/113',
                           'PCB-93/95/98/100/102',
                           'PCB-95/98/102',
                           'Pesticide mix, unspecified',
                           'PBDE 156/169',
                           'BDE-28/33',
                           'meta & para Xylene mix',
                           'meta & para Xylene mix***retired***use m,p-Xylene',
                           'm-cresol mixt. with p-phenol',
                           'Monomeric aluminum',
                           'Oil and Grease',
                           'Oil and grease',
                           'Oil and Grease surface slick/sheen (Y/N) (choice list)',
                           'Oil and Grease, surface slick/sheen - severity (choice list)',
                           'Oil range organics',
                           'Oil range organics C28-C35',
                           'Omeprazole/Esomeprazole mix',
                           'Organic carbon',
                           'Organic compounds',
                           'Organic halides, total -- SWDA NPDWR',
                           'Organic monomeric aluminum',
                           'Organics volatile mix, unspecified',
                           'ortho & para Xylene mix',
                           
                           'Cylindrospermopsins',
                           'Chlordane, technical, and/or chlordane metabolites',
                           'Chlorinated dioxins and furans -- 2,3,7,8 congeners only',
                           '2,4,5-T + Silvex',
                           'Benzene, toluene, ethyl benzene, xylenes mix',
                           'C10-C12 Aliphatics',
                           'C10-C12 Aromatics',
                           'C11-C22 Aromatics',
                           'C19-C36 Aliphatics',
                           'C1-Benzo[a]anthracenes/chrysenes',
                           'C1-C3 Fluorenes',
                           'C1-C4 Chrysenes',
                           'C1-C4 Fluoranthenes',
                           'C1-C4 Phenanthrenes',
                           'C1-Fluoranthenes/pyrenes',
                           'C1-Phenanthrenes/anthracenes',
                           'C23-C32 Hydrocarbons',
                           'C24-C36 Hydrocarbons',
                           'C28 Hydrocarbons',
                           'C2-Benzo[a]anthracenes/chrysenes',
                           'C2-Chrysenes',
                           'C2-Dibenzothiophenes',
                           'C2-Fluoranthenes/pyrenes',
                           'C2-Fluorenes',
                           'C2-Naphthalenes',
                           'C2-Phenanthrenes/anthracenes',
                           'C3-Benzo[a]anthracenes/chrysenes',
                           'C3-Chrysenes',
                           'C3-Dibenzothiophenes',
                           'C3-Fluoranthenes/pyrenes',
                           'C3-Fluorenes',
                           'C3-Naphthalenes',
                           'C3-Phenanthrenes/anthracenes',
                           'C4-Benzo[a]anthracenes/chrysenes',
                           'C4-C6 Aliphatics',
                           'C4-Chrysenes',
                           'C4-Naphthalenes',
                           'C4-Phenanthrenes/anthracenes',
                           'C5-C8 Aliphatics',
                           'C6-C8 Aliphatics',
                           'C8-Alkylphenol',
                           'C8-C10 Aliphatics',
                           'C8-C10 Aromatics',
                           'C9-C10 Aromatics',
                           'C9-C12 Aliphatics',
                           'C9-C18 Aliphatics',
                           "Xylenes mix of m + o + p",
                           "Triazines mixture, unspecified",
                           'C23-C32 Hydrocarbons',"Tannin and Lignin",
                           'Xylenes, m- & p- Mix***retired***use m,p-Xylene',
                           'o-Xylene, mixt. with m-xylene and p-xylene',
                           'Benzene, toluene, ethyl benzene, xylenes mix',
                           'Oil range organics','C1-C3 Fluorenes','C1-C4 Phenanthrenes',
                           'Extended diesel range organics C10-C36',
                           'UV 254','o-Xylene, mixt. with m-xylene and p-xylene',
                           'BHC, .beta.-BHC & .gamma.-BHC mix, unspecified',
                           'Triazines mixture, unspecified',
                           'Pyrethrins','Azinphos-methyl oxygen analog','Microcystin',
                           'Microcystin',
                           'MBAS', 
                           "Hydroxyfluometuron",
                           "Hydroxyphthalazinone",
                           'Microcystin HiLR',
                           'Microcystin HtYR',
                           'Microcystin WR',
                           'BHC, .beta.-BHC & .gamma.-BHC mix, unspecified',
                           'Benzene, 1,2,4,5-tetrachloro- and/or 1,2,3,5-Tetrachlorobenzene',
                           'Extended diesel range organics C10-C36',
                           'Fuel, Diesel Range (C28-C36)',
                           'Diuron Metabolite',
                           'Chrysene + Triphenylene',
                           'Herbicide mix, unspecified',
                           'Inorganic carbon',
                           #metals
                           'Monomeric aluminum','Organic monomeric aluminum',
                           'Stannanetriylium, butyl',
                           'Aluminum, Organic Monomeric (reactive aluminum)',
                           'Aluminum, Organic + Inorganic Monomeric (reactive aluminum)',
                           
                           #ambiguous structure
                           'Octachlorobiphenyl',
                           'Nonachlorobiphenyl',
                           'Dichloroethylene',

                           'Methylfluorene',
                           'Methyldibenzothiophene',
                           'Methylchrysene',
                           'Oil and Grease, surface slick/sheen - severity (choice list)',
                           'Oil and Grease surface slick/sheen (Y/N) (choice list)',
                           'Chlorthal monoacid and diacid degradates','Dissolved microcystins plus nodularins',
                           #radiolabeled
                           'Perylene-d12',
                           'Pyrene-d10',
                           '2,4-D-13C6',
                           'Caffeine-trimethyl-13C3',
                           'Ibuprofen-13C3',
                           'Malathion-D10',
                           'p-Terphenyl-d14',
                           'Phenol-d5',
                           'Atenolol-D7',
                           'Atrazine-d5 (ethyl-d5)',
                           'Alprazolam-D5',
                           'Bromoxynil-13C6',
                           'Caffeine-trimethyl-13C3',
                           'Carbofuran-D3',
                           'Cocaethylene-D3',
                           'Cocaine-D3',
                           'Codeine-d6',
                           'p-Terphenyl-d14',
                           '2,4-D-13C6',
                           '2-Butanone-d5',
                           '2-Chlorophenol-d4',
                           'Pyrene-d10',
                           'Perylene-d12',
                           'Sulfamethoxazole-13C6','Phenol-d5',
                           'Phenobarbital-D5',
                           'Norfluoxetine-D6',
                           'Lorazepam-D4',
                           'Malathion-D10',
                           'Ibuprofen-13C3',
                           'Hydrocodone-D6',
                           'Fluoxetine-D6',
                           'DL-Amphetamine-D5',
                           'Ecgonine Methyl Ester-D3',
                           'Diazepam-D5',
                           'Diazinon-D10',
                           'Benzoylecgonine-D3',
                           'Gemfibrozil-D6',
                           'Diphenhydramine-D3',

                            #! maybe can be registered
                           'Trichloropropene sulfonic acid ethyl ester',
                           'Verrucarin A',
                           'Testosterone sulfate',
                           'Tefluthrin acid pentafluorobenzyl ester',
                           'Tefluthrin acid benzyl ester',
                           'sec-Alachlor sulfonic acid',
                           'Phenol, 2-(1,1-dimethylethyl)-4-methoxy',
                           'O-Ethyl S-methyl S-propyl phosphorodithioate',
                           'Metolachlor hydroxy morpholinone',
                           'Oxymorphone glucuronide',
                           '1,3-Dinitroso-5-nitro-1,3,5-triazacyclohexane',
                           '17-beta-Estradiol 17-sulfate',
                           '1-Nitroso-3,5-dinitro-1,3,5-triazacyclohexane',
                           'Hydroxytebuthiuron',
                           'Geosmin',
                           'Diclofenac, 4-hydroxy'
                           ]
#! not yet checked
noncurated_CharacteristicNames = ['2-(1-Hydroxyethyl)-6-methylaniline',
                                  '2-(4-Octylphenoxy)ethanol',
                                  '2-(Methylsulfonyl)-4-(trifluoromethyl)benzoic acid (diketonitrile-isoxaflutole benzoic acid analog)',
                                  '2-(p-tert-Butylphenoxy)cyclohexanol',
                                  "2,2',3,3',4',5,6-Heptachlorobiphenyl",
                                  "2,2',3,4,4',5,5',6-OCBDE",
                                  "2,2',4,4',5',6-HXBDE",
                                  "2,2'-[1,2 ethanediylbis(oxy)]bis-ethanol diacetate",
                                  '2,3,3-Trichloro-2-propene-1-sulfonic acid (sodium salt)',
                                  "2,3',4',6-TEBDE",
                                  '2,3,6-Trichloro-5-cyano-4-hydroxybenzamide',
                                  '2,3,7,8-Tetrachlorodibenzo-p-dioxin, TEQ',
                                  '2,4,5,6-Tetrachloroisophthalamide',
                                  '2,4,6-TRBDE',
                                  '2,4-DIBDE',
                                  '2,6-DIBDE',
                                  '2-[(2-Ethyl-6-methylphenyl)amino]-2-oxoethanesulfonic acid',
                                  '2-Isobutyl-3-methoxypyraxine',
                                  '2-Methylisoborneol',
                                  '2-MOBDE',
                                  '3-(Trifluoromethyl)phenylurea',
                                  "3,3',4,4',5-PEBDE",
                                  "3,3',4,4'-TEBDE",
                                  '3-MOBDE',
                                  '4(3H)-Quinazolinone, 2-methyl-3-(2-methylphenyl)',
                                  '4-Chlorophenyl ether',
                                  '4-Epianhydrochlortetracycline',
                                  '4H-1-Benzopyran-4-one, 5,7-dihydroxy-3-(4-hydroxyphenyl)',
                                  '4H-1-Benzopyran-4-one, 5,7-dihydroxy-3-(4-methoxyphenyl)',
                                  '4-Hydroxy molinate',
                                  '4-Hydroxyhexazinone A',
                                  '7-Hydroxyquetiapine',
                                  'AB-PINACA 4-hydroxypentyl',
                                  'Acetamide, 2-(diethylamino)-N-(2,6-dimethylphenyl)',
                                  'alpha.,.alpha.-Dimethylphenethylamine',
                                  'alpha.-Amino-2,3-dihydro-5-methyl-3-oxo-4-isoxazolepropanoic acid',
                                  'alpha.-Apo-oxytetracycline',
                                  'Amitriptyline (+/-)-E-10-hydroxylated',
                                  'Androsterone glucuronide',
                                  'Androsterone sulfate',
                                  'Anhydroerthromycin',
                                  'Anhydrotetracycline',
                                  "Benzene,1,1'-(chloroethenylidene)bis(4-chloro-)",
                                  'Benzeneacetic acid, .alpha.-methyl-4-(2-methylpropyl)',
                                  'Benzeneethanamine, .alpha.-methyl-, (.alpha.S)-, sulfate (2:1)',
                                  'Benzenemethanol, .alpha.-[(1S)-1-(methylamino)ethyl]-, (.alpha.S)',
                                  'Benzenesulfonamide, 4-amino-N-2-pyrimidinyl',
                                  'Benzo[b,k]fluoranthene',
                                  'Benzoylecgonine hydrate',
                                  'beta.-Apo-oxytetracycline',
                                  'Bicyclopyrone SYN503870',
                                  'Bis(hexachlorocyclopentadieno) cyclooctane',
                                  'Brombuterol',
                                  'Bromoconazole',
                                  'Bromodichloropropane',
                                  'Butachlor ESA',
                                  'Carbamazepine 10,11 epoxide',
                                  'Carboxy molinate',
                                  'Chlorosulfonamide acid',
                                  'cis-Cyhalothric acid',
                                  'Cyanazine acid',
                                  'Cyclopropanecarboxylic acid, 3-(2,2-dichloroethenyl)-2,2-dimethyl-, methyl ester, (1R,3R)-rel',
                                  'D9 (+/-)11-nor-9-carboxy-delta-THC',
                                  "DDE, o,o'-, o,p'-, and p,p'- isomers",
                                  "DDE, o,p'- and p,p'- isomers",
                                  "DDT, o,p'- and p,p'- isomers",
                                  'Dechloroacetochlor',
                                  'Dechloroalachlor',
                                  'Dechlorodimethenamid',
                                  'Dechlorofipronil',
                                  'Dechlorometolachlor',
                                  'Deethylcyanazine',
                                  'Deethylcyanazine acid',
                                  'Deethylcyanazine amide',
                                  'Dehydroaripiprazole hydrochloride',
                                  'Deiodo flubendiamide',
                                  'Deisopropyl prometryn',
                                  'Deisopropylprometryn',
                                  'Demethyl hexazinone B',
                                  'Demethylfluometuron',
                                  'Desmethylvenlafaxine',
                                  'Desulfinylfipronil amide',
                                  'Di(2-ethoxylhexyl) phthalate',
                                  'Dibutylstannanediylium',
                                  'Didemethyl hexazinone F',
                                  'Didemethyl tebuthiuron',
                                  'Diethylstilbestrol glucuronide (DES-GlcA)',
                                  'Dimethenamid sulfinylacetic acid',
                                  'Disulfoton oxon sulfone',
                                  'Disulfoton oxon sulfoxide',
                                  'Epi-chlorotetracycline',
                                  'Epi-iso-chlorotetracycline',
                                  'Epi-oxytetracycline',
                                  'EPTC degradate R248722',
                                  'Equilenin sulfate',
                                  'Erythromycin-anhydro',
                                  'Estriol 17-sulfate',
                                  'Estriol 3-glucuronide',
                                  'Estriol 3-sulfate',
                                  'Estrone 3-sulfate',
                                  'Estrone glucuronide',
                                  'Ethanamine, 2-(diphenylmethoxy)-N,N-dimethyl',
                                  'Ethane, 1-chloro-1,1-difluoro',
                                  'Ethanol, 2-(4-nonylphenoxy)',
                                  'Ethanone, 1-(2,3-dihydro-1,1,2,3,3,6-hexamethyl-1H-inden-5-yl)',
                                  'Ethanone, 1-[6-(1,1-dimethylethyl)-2,3-dihydro-1,1-dimethyl-1H-inden-4-yl]',
                                  'Ethenylestradiol 3-glucuronide',
                                  'Ethinylestradiol 3-sulfate',
                                  'Fipronil sulfonate',
                                  'Hexachlorobutene',
                                  'Hydromorphone hydrochloride',
                                  'Hydroxy didemethyl fluometuron',
                                  'Hydroxy monodemethyl fluometuron',
                                  'Hydroxyacetochlor',
                                  'Hydroxybupropion',
                                  'Hydroxydiazinon',
                                  'Hydroxydimethenamid']
    
no_SMILES = ['Octylphenol',
'Nonylphenol',
'EDDP perchlorate',
'Aroclor 1016',
'Aroclor 1221',
'Aroclor 1232',
'Aroclor 1242',
'Aroclor 1248',
'Aroclor 1254',
'Aroclor 1260',
'Aroclor 1262',
'Aroclor 1268',
'Asbestos',
'Butylated hydroxyanisole',
'Oxytetracycline calcium',
'Chlorophenol',
'Chlorotoluene',
'Cresol',
'Demeton',
'Demeton-methyl',
'Dichlorobenzene',
'Cacodylic acid',
'Dimethylnaphthalene',
'Dinitrotoluene',
'Emamectin benzoate',
'1,2-Dichloroethane-d4',
'Fenbutatin-oxide',
'Silvex isooctyl ester',
'Fentin',
'Heptachlorodibenzo-p-dioxin',
'Hexabromodiphenyl ether',
'Hexachlorocyclohexane',
'Hexachlorodibenzo-p-dioxin',
'Hexachlorodibenzofuran',
'Hydrocarbons',
'C12 Hydrocarbons',
'Imazamethabenz-methyl',
'Kerosene',
'm,p-Xylene',
'Methylmercury(1+)',
'Monomethylarsonate',
'Bisphenol F',
'Methylnaphthalene',
'Petroleum spirits',
'Nonylphenol',
'Nonabromophenoxybenzene',
'Petroleum phenols',
'FireMaster BP 6',
'Polybrominated biphenyls',
'Polychlorinated biphenyls',
'Chlorinated naphthalenes',
'Pyrethrins',
'Quetiapine fumarate',
'Chlordane, technical',
'Terpineol',
'Tetrabutyltin',
'Tetrachlorodibenzo-p-dioxin',
'Tetrachloroethane',
'Tetrachlorophenol',
'Tolyl triazole',
'Toxaphene',
'Trichlorobenzene',
'Trichlorophenol',
'Trihalomethanes',
'Virginiamycin',
'Xylene',
'Zinc phosphide',
'Ziram']                                 
                                  
high_freq_chem_names = ['Atrazine',
'Metolachlor',
'Alachlor',
'Simazine',
'Acetochlor',
'Metribuzin',
'Cyanazine',
'Prometon',
'2-Chloro-4-isopropylamino-6-amino-s-triazine',
'Tetrachloroethylene',
'Diazinon',
'Chlorpyrifos',
'Carbaryl',
'Trichloroethylene',
'Carbofuran',
'Trifluralin',
'Malathion',
'Benzene',
'Toluene',
'Vinyl chloride',
'1,1-Dichloroethylene',
'Ethylbenzene',
'p-Dichlorobenzene',
'S-Ethyl dipropylthiocarbamate',
'Chloroform',
'trans-1,2-Dichloroethylene',
'Pendimethalin',
'Methylene chloride',
'Carbon tetrachloride',
'cis-1,2-Dichloroethylene']

# A. load data
# A.1. load and filter station data



"https://www.waterqualitydata.us/data/Station/search?countrycode=US&characteristicType=Organics%2C%20Other&characteristicType=Organics%2C%20PCBs&characteristicType=Organics%2C%20Pesticide&characteristicType=PFAS%2CPerfluorinated%20Alkyl%20Substance&characteristicType=PFOA%2C%20Perfluorooctanoic%20Acid&startDateLo=01-01-2018&startDateHi=02-01-2018&mimeType=csv"

url_prefix = "https://www.waterqualitydata.us/data/"

url_part_2 = "/search?countrycode=US&sampleMedia=water&sampleMedia=Water&characteristicType=Organics%2C%20Other&characteristicType=Organics%2C%20PCBs&characteristicType=Organics%2C%20Pesticide&characteristicType=PFAS%2CPerfluorinated%20Alkyl%20Substance&characteristicType=PFOA%2C%20Perfluorooctanoic%20Acid&startDateLo="

url_post_date_01 = "-01-"

url_post_date_02 = "&startDateHi="

url_post_date_03 = "-31-"

url_end = "&mimeType=csv&zip=yes"

dates = np.arange(2008,2019,1) #!2008
months = [("01","03"),("04","07"),("08","10"),("11","12")]
seasons = [("10","11","12"),
           ("01","02","03"),
           ('04','05','06'),
           ('07','08','09')]

bad_station_types = ['Land Flood Plain',
                     'Spigot / Faucet',
                     'Facility: Laboratory or sample-preparation area',
                     'Wetland Undifferentiated',
                     'Spigot / Faucet',
                     'BEACH Program Site-Estuary',
                     'Mine/Mine Discharge',
                     'Land: Outcrop',
                     'Subsurface: Groundwater drain',
                     'Facility: Wastewater land application',
                     'Mine Pit',
                     'Subsurface: Cave',
                     'Cave',
                     'Atmosphere',
                     'Local Air Monitoring Station',
                     'Landfill',
                     'BEACH Program Site-Ocean',
                     'Gallery',
                     'Estuary',
                     'Land',
                     'Facility: Cistern',
                     'Facility: Waste injection well',
                     'Facility Industrial',
                     'Facility Public Water Supply (PWS)',
                     'Ocean: Coastal',
                     'Subsurface: Tunnel, shaft, or mine',
                     'Stream: Tidal stream',
                     'Waste Pit',
                     'Well',
                     'Facility: Septic system',
                     'Leachate-SamplePoint',
                     'Well: Multiple wells',
                     'Pipe, Unspecified Source',
                     'Facility: Wastewater sewer',
                     'Facility: Combined sewer',
                     'Waste Sewer',
                     'Facility Municipal Sewage (POTW)',
                     'Borehole',
                     'Ocean',
                     "Borehole"]

#! did not do this yet - left them in to compare values
unique_MonitoringLocationTypeName = ['Land Flood Plain',
 'Facility: Laboratory or sample-preparation area',
 'Wetland Undifferentiated',
 'Spigot / Faucet',
 'Great Lake',
 'BEACH Program Site-Estuary',
 'Mine/Mine Discharge',
 'Well: Hyporheic-zone well',
 'Land: Outcrop',
 'Subsurface: Groundwater drain',
 'Wetland Lacustrine-Emergent',
 'Other-Surface Water',
 'Storm Sewer',
 'Facility: Wastewater land application',
 'River/Stream',
 'Mine Pit',
 'Ocean',
 'Facility Other',
 'Aggregate surface-water-use',
 'Subsurface: Cave',
 'Facility: Pavement',
 'Aggregate groundwater use',
 'Atmosphere',
 'Canal Transport',
 'Facility: Storm sewer',
 'Stream',
 'Landfill',
 'Land',
 'Facility: Cistern',
 'BEACH Program Site-River/Stream',
 'Wetland Estuarine-Emergent',
 'Spring',
 'Facility: Diversion',
 'Facility: Combined sewer',
 'BEACH Program Site-Lake',
 'Combined Sewer',
 'Ocean: Coastal',
 'Land Runoff',
 'Well: Collector or Ranney type well',
 'Wetland',
 'Facility: Wastewater sewer',
 'Borehole',
 'Riverine Impoundment',
 'River/Stream Intermittent',
 'Subsurface',
 'Land: Excavation',
 'Facility: Water-use establishment',
 'Wetland Palustrine-Forested',
 'Waste Sewer',
 'River/stream Effluent-Dominated',
 'Pond-Stormwater',
 'Channelized Stream',
 'Facility: Field, Pasture, Orchard, or Nursery',
 'CERCLA Superfund Site',
 'Stream: Canal',
 'Well: Test hole not completed as a well',
 'Reservoir',
 'Facility: Waste injection well',
 'Lake',
 'Other-Ground Water',
 'Local Air Monitoring Station',
 'Wetland Palustrine Pond',
 'Stream: Ditch',
 'Facility Public Water Supply (PWS)',
 'Wetland Riverine-Emergent',
 'River/Stream Perennial',
 'Lake, Reservoir, Impoundment',
 'Cave',
 'Facility Municipal Sewage (POTW)',
 'Seep',
 'Facility Privately Owned Non-industrial',
 'Facility: Water-distribution system',
 'Constructed Wetland',
 'River/Stream Ephemeral',
 'Land: Sinkhole',
 'Pond-Stock',
 'Subsurface: Tunnel, shaft, or mine',
 'Stream: Tidal stream',
 'Waste Pit',
 'Canal Drainage',
 'Well',
 'Facility: Septic system',
 'Leachate-SamplePoint',
 'Well: Multiple wells',
 'Pipe, Unspecified Source',
 'Wetland Palustrine-Shrub-Scrub',
 'BEACH Program Site-Ocean',
 'Gallery',
 'Wetland Palustrine-Emergent',
 'Estuary',
 'Subsurface: Unsaturated zone',
 'Canal Irrigation',
 'Facility Industrial',
 'Facility: Outfall']

bad_fips = [2,3,7,14,15,43,
            52,60,64,66,67,68,69,70,71,72,73,74,76,78,
            81,84,86,89]
#stations = stations[~stations.MonitoringLocationTypeName.isin(bad_station_types)]

#! stratify on MonitoringLocationTypeName

#too_big = []
#data_len = []

station_df = pd.DataFrame(columns=['OrganizationIdentifier',
                                 'OrganizationFormalName',
                                 'MonitoringLocationIdentifier',
                                 'MonitoringLocationName',
                                 'MonitoringLocationTypeName',
                                 'MonitoringLocationDescriptionText',
                                 'HUCEightDigitCode',
                                 'DrainageAreaMeasure/MeasureValue',
                                 'DrainageAreaMeasure/MeasureUnitCode',
                                 'ContributingDrainageAreaMeasure/MeasureValue',
                                 'ContributingDrainageAreaMeasure/MeasureUnitCode',
                                 'LatitudeMeasure',
                                 'LongitudeMeasure',
                                 'SourceMapScaleNumeric',
                                 'HorizontalAccuracyMeasure/MeasureValue',
                                 'HorizontalAccuracyMeasure/MeasureUnitCode',
                                 'HorizontalCollectionMethodName',
                                 'HorizontalCoordinateReferenceSystemDatumName',
                                 'VerticalMeasure/MeasureValue',
                                 'VerticalMeasure/MeasureUnitCode',
                                 'VerticalAccuracyMeasure/MeasureValue',
                                 'VerticalAccuracyMeasure/MeasureUnitCode',
                                 'VerticalCollectionMethodName',
                                 'VerticalCoordinateReferenceSystemDatumName',
                                 'CountryCode',
                                 'StateCode',
                                 'CountyCode',
                                 'AquiferName',
                                 'FormationTypeText',
                                 'AquiferTypeName',
                                 'ConstructionDateText',
                                 'WellDepthMeasure/MeasureValue',
                                 'WellDepthMeasure/MeasureUnitCode',
                                 'WellHoleDepthMeasure/MeasureValue',
                                 'WellHoleDepthMeasure/MeasureUnitCode',
                                 'ProviderName']
                        )

problems = []    
for d in dates:
    for m in months:
        this_url = url_prefix + "Station" + url_part_2 + str(m[0]) + url_post_date_01 + str(d) + url_post_date_02 + str(m[1]) + url_post_date_03 + str(d) + url_end
        try:
            data = pd.read_csv(this_url, compression='zip', sep=',', low_memory=False)    
            station_df = pd.concat([station_df, data])
            station_df.drop_duplicates(inplace=True)
            station_df = station_df[~station_df.MonitoringLocationTypeName.isin(bad_station_types)]
            station_df = station_df[~station_df.StateCode.isin(bad_fips)]
            print("Loading " + str(m[0]) + " to " + str(m[1]) + " " + str(d))
        except: 
            problems.append((m,d))
            print(str(m[0]) + " to " + str(m[1]) + " " + str(d) + " broke :(")

#! find a few HUCs
station_info_df = station_df[['MonitoringLocationIdentifier',
                              'MonitoringLocationTypeName',
                              'HUCEightDigitCode']]

# station_df.to_pickle("station_df.pkl") save for later analysis

# A.2. load and filter activity data, result data, flag bloq
#! may want to remove flow-through?
# strat on depth or altitude?

all_activity_columns = ['OrganizationIdentifier',
                        'OrganizationFormalName',
                        'ActivityIdentifier',
                        'ActivityTypeCode',
                        'ActivityMediaName',
                        'ActivityMediaSubdivisionName',
                        'ActivityStartDate',
                        'ActivityStartTime/Time',
                        'ActivityStartTime/TimeZoneCode',
                        'ActivityEndDate',
                        'ActivityEndTime/Time',
                        'ActivityEndTime/TimeZoneCode',
                        'ActivityRelativeDepthName',
                        'ActivityDepthHeightMeasure/MeasureValue',
                        'ActivityDepthHeightMeasure/MeasureUnitCode',
                        'ActivityDepthAltitudeReferencePointText',
                        'ActivityTopDepthHeightMeasure/MeasureValue',
                        'ActivityTopDepthHeightMeasure/MeasureUnitCode',
                        'ActivityBottomDepthHeightMeasure/MeasureValue',
                        'ActivityBottomDepthHeightMeasure/MeasureUnitCode',
                        'ProjectIdentifier',
                        'ActivityConductingOrganizationText',
                        'MonitoringLocationIdentifier',
                        'ActivityCommentText',
                        'SampleAquifer',
                        'HydrologicCondition',
                        'HydrologicEvent',
                        'ActivityLocation/LatitudeMeasure',
                        'ActivityLocation/LongitudeMeasure',
                        'ActivityLocation/SourceMapScaleNumeric',
                        'ActivityLocation/HorizontalAccuracyMeasure/MeasureValue',
                        'ActivityLocation/HorizontalAccuracyMeasure/MeasureUnitCode',
                        'ActivityLocation/HorizontalCollectionMethodName',
                        'ActivityLocation/HorizontalCoordinateReferenceSystemDatumName',
                        'AssemblageSampledName',
                        'CollectionDuration/MeasureValue',
                        'CollectionDuration/MeasureUnitCode',
                        'SamplingComponentName',
                        'SamplingComponentPlaceInSeriesNumeric',
                        'ReachLengthMeasure/MeasureValue',
                        'ReachLengthMeasure/MeasureUnitCode',
                        'ReachWidthMeasure/MeasureValue',
                        'ReachWidthMeasure/MeasureUnitCode',
                        'PassCount',
                        'NetTypeName',
                        'NetSurfaceAreaMeasure/MeasureValue',
                        'NetSurfaceAreaMeasure/MeasureUnitCode',
                        'NetMeshSizeMeasure/MeasureValue',
                        'NetMeshSizeMeasure/MeasureUnitCode',
                        'BoatSpeedMeasure/MeasureValue',
                        'BoatSpeedMeasure/MeasureUnitCode',
                        'CurrentSpeedMeasure/MeasureValue',
                        'CurrentSpeedMeasure/MeasureUnitCode',
                        'ToxicityTestType',
                        'SampleCollectionMethod/MethodIdentifier',
                        'SampleCollectionMethod/MethodIdentifierContext',
                        'SampleCollectionMethod/MethodName',
                        'SampleCollectionMethod/MethodQualifierTypeName',
                        'SampleCollectionMethod/MethodDescriptionText',
                        'SampleCollectionEquipmentName',
                        'SampleCollectionMethod/SampleCollectionEquipmentCommentText',
                        'SamplePreparationMethod/MethodIdentifier',
                        'SamplePreparationMethod/MethodIdentifierContext',
                        'SamplePreparationMethod/MethodName',
                        'SamplePreparationMethod/MethodQualifierTypeName',
                        'SamplePreparationMethod/MethodDescriptionText',
                        'SampleContainerTypeName',
                        'SampleContainerColorName',
                        'ChemicalPreservativeUsedName',
                        'ThermalPreservativeUsedName',
                        'SampleTransportStorageDescription',
                        'BinaryObjectFileName',
                        'BinaryObjectFileTypeCode',
                        'ActivityFileUrl',
                        'ActivityMetricUrl',
                        'ActivityGroupUrl',
                        'ProviderName']
 
bad_activity_columns = ['SampleContainerTypeName',
                        'SamplePreparationMethod/MethodName',
                        'SamplePreparationMethod/MethodQualifierTypeName',
                        'ThermalPreservativeUsedName',
                        'SampleContainerColorName',
                        'ChemicalPreservativeUsedName',
                        'SamplePreparationMethod/MethodDescriptionText',
                        'BinaryObjectFileName',
                        'BinaryObjectFileTypeCode',
                        'ActivityFileUrl',
                        'BoatSpeedMeasure/MeasureValue',
                        'ToxicityTestType',
                        'CurrentSpeedMeasure/MeasureUnitCode',
                        'CurrentSpeedMeasure/MeasureValue',
                        'BoatSpeedMeasure/MeasureUnitCode',
                        'NetMeshSizeMeasure/MeasureUnitCode',
                        'NetMeshSizeMeasure/MeasureValue',
                        'NetSurfaceAreaMeasure/MeasureUnitCode',
                        'NetSurfaceAreaMeasure/MeasureValue',
                        'NetTypeName',
                        'PassCount',
                        'ReachLengthMeasure/MeasureValue',
                        'SamplingComponentPlaceInSeriesNumeric',
                        'ReachLengthMeasure/MeasureUnitCode',
                        'ReachWidthMeasure/MeasureValue',
                        'ReachWidthMeasure/MeasureUnitCode',
                        'ActivityMetricUrl',
                        'ActivityGroupUrl']

bad_activity_types = ['Quality Control Sample-Equipment Blank',                   
                      'Quality Control Sample-Field Blank',
                      'Quality Control Sample-Trip Blank',
                      'Sample-Negative Control',
                      'Quality Control Sample-Field Spike',
                      'Quality Control Sample-Field Ambient Conditions Blank',
                      'Quality Control Sample-Lab Matrix Spike',
                      'Quality Control Sample-Lab Matrix Spike Duplicate',
                      'Quality Control Sample-Lab Control Sample/Blank Spike',
                      'Quality Control Sample-Lab Spike',
                      'Quality Control Field Sample Equipment Rinsate Blank',
                      'Quality Control Sample-Spike Solution']

bad_ActivityMediaSubdivisionNames = ['Ocean Water',
                                     #'Ground Water',
                                     #'Groundwater',
                                     #'Mixing Zone, Zone of Initial Dilution',
                                     #'Elutriation', 
                                     #'Wet Fall Material',
                                     #'Wastewater Treatment Plant Influent',
                                     'Finished Water'] #! have a chemist check a list of all activity types

some_ActivityMediaSubdivisionNames = ['Surface Water', 'Ground Water', 'Groundwater', 'Estuary',
       'Water', 'Septic Effluent', 'Ocean Water', 'Leachate',
       'Industrial Waste', 'Wastewater Treatment Plant Effluent',
       'Bulk deposition', 'Wet Fall Material', 'Mixing Zone',
       'Stormwater', 'Mixing Zone, Zone of Initial Dilution',
       'Landfill effluent', 'Industrial Effluent', 'Interstitial',
       'Elutriation', 'Wastewater Treatment Plant Influent',
       'Finished Water']

bad_water_types = ['Bulk deposition', #! not currently used. what is this field?
                   'Deionized Water',
                   'Elutriation',
                   'Finished Water',
                   'Ground Water',
                   'Groundwater',
                   'Hyporheic zone',
                   'Industrial Waste',
                   'Interstitial',
                   'Interstitial Water',
                   'Municipal Waste',
                   'Ocean Water',
                   'Pore water',
                   'Snowmelt',
                   'Surface Water Sediment',
                   'Wet Fall Material']

bad_result_columns = ['TaxonomistAccreditationIndicator','DataLoggerLine',
           'TaxonomistAccreditationAuthorityName',
           'CellShapeName', 'HabitName', 'VoltismName',
       'TaxonomicPollutionTolerance', 'TaxonomicPollutionToleranceScaleText',
       'TrophicLevelName', 'FunctionalFeedingGroupName',
       'TaxonomicDetailsCitation/ResourceTitleName',
       'FrequencyClassInformationUrl',
       'TaxonomicDetailsCitation/ResourceCreatorName',
       'TaxonomicDetailsCitation/ResourceSubjectText',
       'TaxonomicDetailsCitation/ResourcePublisherName',
       'TaxonomicDetailsCitation/ResourceDate','MethodSpecificationName',
       'TaxonomicDetailsCitation/ResourceIdentifier',
       'BiologicalIntentName', 'BiologicalIndividualIdentifier',
       'SubjectTaxonomicName', 'UnidentifiedSpeciesIdentifier',
       'SampleTissueAnatomyName', 'GroupSummaryCountWeight/MeasureValue',
       'GroupSummaryCountWeight/MeasureUnitCode', 'CellFormName',
                      'ResultSamplingPointName',
                      'ResultTimeBasisText',
                      'ResultTemperatureBasisText',
                      'LaboratoryAccreditationIndicator',
                      'LaboratoryAccreditationAuthorityName']

# harmonized ResultDetectionConditionText as binary BLOQ
# Systematic Contamination sounds weird, but spot check looks okay

#! check BLOQ labeling -- some seem to have gotten through
bloq_rdct = ['Present Below Quantification Limit',
        'Below Method Detection Limit',
        'Below Reporting Limit',
        'Detected Not Quantified',
        '*Non-detect',
        'Not Detected']

bloq_rmv = ['-',
            '?',
            '*',
            '<',
            'BDL',
            'n.d.',
            'nd',
            'ND',
            'No Data',
            'non detect',
            'Non-detect',
            'Not Detected',
            'Q',
            0            ] #! check No Data, Q
# checked ResultLaboratoryCommentText contains below, didn't add value -- can look for other data in this field though
bloq_sbc = 'Maximum' #'Minimum'? StatisticalBaseCode
#! ResultLaboratoryCommentText contains 'value extrapolated at low' should be added

def is_bloq(row):
    if row['ResultDetectionConditionText'] in bloq_rdct:
        return 1
    elif str(row['ResultMeasureValue']).startswith(('-',
            '?',
            '*',
            '<',
            'BDL',
            'n',
            'N',
            'Q')):
        return 1
    elif row['ResultMeasureValue'] == 0:
        return 1
    elif row['StatisticalBaseCode'] == bloq_sbc :
        return 1
    else: return 0

technical_min = ['Method Detection Level',
 'Elevated Detection Limit',
 'Instrument Detection Level',
 'Estimated Detection Level',
 'Method Detection Limit (MDL)']

technical_qnt = ['Practical Quantitation Limit','Lower Quantitation Limit']

reporting_min = ['Minimum Reporting Level',
 'Lower Reporting Limit',
 'Laboratory Reporting Level',
 'Historical Lower Reporting Limit']

#! strip lt sign from loq values
#! check ones where loq = res
def calc_loq_val(row):
    try:
        this_loq_val = float(row['DetectionQuantitationLimitMeasure/MeasureValue'])
        if row['DetectionQuantitationLimitMeasure/MeasureUnitCode'] in ['ug/l','ug/L','ppb']:
            return this_loq_val
        elif row['DetectionQuantitationLimitMeasure/MeasureUnitCode'] in ['ng/l','ng/L','ppt']:
            return this_loq_val/1000
        elif row['DetectionQuantitationLimitMeasure/MeasureUnitCode'] in ['mg/l','mg/L']:
            return this_loq_val*1000
        else: return float('NaN')   
    except: float('NaN') 
        
def calc_res_val(row):
    try: 
        if row['bloq'] == 0:
            if row['ResultMeasure/MeasureUnitCode'] in ['ug/l','ug/L','ppb']:
                return float(row['ResultMeasureValue'])
            elif row['ResultMeasure/MeasureUnitCode'] in ['ng/l','ng/L','ppt']:
                return float(row['ResultMeasureValue'])/1000
            elif row['ResultMeasure/MeasureUnitCode'] in ['mg/l','mg/L']:
                return float(row['ResultMeasureValue'])*1000
        else: return float('NaN') 
    except: return float('NaN') 
#! make new calc method function

bad_units = ['ng/POCIS','ng','ug']

# phase codes: 1: suspended, 2: dissolved, 3: bulk
# by USGSPCode
usgs_codes = pd.read_excel("USGS_parameter_codes.xlsx") # https://nwis.waterdata.usgs.gov/usa/nwis/pmcodes

ss_codes = usgs_codes[usgs_codes['Parameter Name/Description'].str.contains('suspended sediment') == True]['Parameter Code']
dissolved_codes = usgs_codes[usgs_codes['Parameter Name/Description'].str.contains(' filtered') == True]['Parameter Code']
bulk_codes = usgs_codes[usgs_codes['Parameter Name/Description'].str.contains('unfiltered') == True]['Parameter Code']

#water_lst = ['water','suspended sediment','solids']
#water_terms = '|'.join(water_lst)
#water_codes = usgs_codes[usgs_codes['Parameter Name/Description'].str.contains(water_terms) == True]['Parameter Code']

#len(result[(result.method_1 == 3) & (result.method_2 == 3)])
#len(result[(result.method_1 == 2) & (result.method_2 == 2)])
#len(result[(result.method_1 == 2) & (result.method_2 == 3)])
#len(result[(result.method_1 == 3) & (result.method_2 == 2)])

#methods codes to exclude
air_codes = usgs_codes[usgs_codes['Parameter Name/Description'].str.contains(', air,') == True]['Parameter Code']
bed_sed_codes = usgs_codes[usgs_codes['Parameter Name/Description'].str.contains(', bed sediment,') == True]['Parameter Code']
biota_codes = usgs_codes[usgs_codes['Parameter Name/Description'].str.contains(', biota, ') == True]['Parameter Code']
solid_codes = usgs_codes[usgs_codes['Parameter Name/Description'].str.contains(', solids, ') == True]['Parameter Code']

#! some only have a filtered method, so maybe they are always dissolved?

#! plot the compound for which there are the most frac_text
# by ResultSampleFractionText (left out  'Extractable', 'Total Residual', 'Recoverable')

# by ResultAnalyticalMethod/MethodIdentifier
# from https://water.usgs.gov/XML/NWIS/5.0/ReferenceLists/MethodList.xml
# removed invalid chars such as &, < (replaced with lt)
root = ET.parse('MethodList.xml').getroot()
ss_lst = []
dissolved_lst = []
bulk_lst = []
none_lst = []

for i in range(len(root)):
    try:
        if root[i][2].text.find('water') > 0:
            if root[i][2].text.find('unfiltered') > 0:
                bulk_lst.append((root[i][0].text,root[i][1].text))
            elif root[i][2].text.find(' filtered') > 0:
                dissolved_lst.append((root[i][0].text,root[i][1].text))
            else:
                continue
        elif root[i][2].text.find('suspended sediment') > 0:
            ss_lst.append((root[i][0].text,root[i][1].text))            
        else: 
            none_lst.append((root[i][0].text,root[i][1].text))
    except:
        continue       

#! also check ResultAnalyticalMethod/MethodUrl & MethodDescriptionText
# identify method: 'ResultAnalyticalMethod/MethodName'
# ResultMeasure/MeasureUnitCode as original
##--##--##--##--##

#bad_ResultWeightBasisText = 'Dry'

#types = ['Station',"Result","Activity","ActivityMetric","ResultDetectionQuantitationLimit"]
#activity_add_on_part_two = "siteType=Facility&"
            
bad_unit_types = ["ng/SPMD", "%"] #'ResultMeasure/MeasureUnitCode'

bad_sample_types = ['Bed Sediment']#'ResultSampleFractionText'            

dtxsids = pd.read_excel('dtxsids_from_dashboard.xls')
chem_lst = dtxsids.QSAR_READY_SMILES.unique()
chem_id_lst = []
for i in range(len(chem_lst)):
    chem_id_lst.append((i,chem_lst[i]))
chem_id_df = pd.DataFrame(chem_id_lst, columns=['chem_id','QSAR_READY_SMILES'])

dtxsids = pd.merge(dtxsids,
                  chem_id_df,
                  how='left',
                  on='QSAR_READY_SMILES')
chem_id_lst = list(dtxsids[['INPUT', 'chem_id']].values)

#! make this faster
def assign_chem_id(row):
    this_chem = row['CharacteristicName_x']
    for l in chem_id_lst:
        if l[0] == this_chem:
            return l[1]
    else: return float('NaN')    

def get_method(row):
    if (row.method_1 == 1 or row.method_2 == 1 or row.method_3 == 1) == True:
        return 1
    elif (row.method_1 == 2 or row.method_2 == 2 or row.method_3 == 2) == True:
        return 2
    elif (row.method_1 == 3 or row.method_2 == 3 or row.method_3 == 3) == True:
        return 3
    else: return 0
                      
problems = [] 
   
activity_df = pd.DataFrame(columns = ['ActivityBottomDepthHeightMeasure/MeasureUnitCode',
 'ActivityBottomDepthHeightMeasure/MeasureValue',
 'ActivityCommentText',
 'ActivityConductingOrganizationText',
 'ActivityDepthAltitudeReferencePointText',
 'ActivityDepthHeightMeasure/MeasureUnitCode',
 'ActivityDepthHeightMeasure/MeasureValue',
 'ActivityEndDate',
 'ActivityEndTime/Time',
 'ActivityEndTime/TimeZoneCode',
 'ActivityIdentifier',
 'ActivityLocation/HorizontalAccuracyMeasure/MeasureUnitCode',
 'ActivityLocation/HorizontalAccuracyMeasure/MeasureValue',
 'ActivityLocation/HorizontalCollectionMethodName',
 'ActivityLocation/HorizontalCoordinateReferenceSystemDatumName',
 'ActivityLocation/LatitudeMeasure',
 'ActivityLocation/LongitudeMeasure',
 'ActivityLocation/SourceMapScaleNumeric',
 'ActivityMediaName',
 'ActivityMediaSubdivisionName',
 'ActivityRelativeDepthName',
 'ActivityStartDate',
 'ActivityStartTime/Time',
 'ActivityStartTime/TimeZoneCode',
 'ActivityTopDepthHeightMeasure/MeasureUnitCode',
 'ActivityTopDepthHeightMeasure/MeasureValue',
 'ActivityTypeCode',
 'AssemblageSampledName',
 'CollectionDuration/MeasureUnitCode',
 'CollectionDuration/MeasureValue',
 'HydrologicCondition',
 'HydrologicEvent',
 'MonitoringLocationIdentifier',
 'OrganizationFormalName',
 'OrganizationIdentifier',
 'ProjectIdentifier',
 'ProviderName',
 'SampleAquifer',
 'SampleCollectionEquipmentName',
 'SampleCollectionMethod/MethodDescriptionText',
 'SampleCollectionMethod/MethodIdentifier',
 'SampleCollectionMethod/MethodIdentifierContext',
 'SampleCollectionMethod/MethodName',
 'SampleCollectionMethod/MethodQualifierTypeName',
 'SampleCollectionMethod/SampleCollectionEquipmentCommentText',
 'SamplePreparationMethod/MethodIdentifier',
 'SamplePreparationMethod/MethodIdentifierContext',
 'SampleTransportStorageDescription',
 'SamplingComponentName'])

all_CharacteristicNames = []

for d in dates:
    for m in months:
        activity_url = url_prefix + 'Activity' + url_part_2 + str(m[0]) + url_post_date_01 + str(d) + url_post_date_02 + str(m[1]) + url_post_date_03 + str(d) + url_end + "&dataProfile=activityAll"
        try:
            activities = pd.read_csv(activity_url, compression='zip', sep=',', low_memory=False)
            # remove unused columns
            activities = activities.drop(bad_activity_columns, axis=1)
            # remove unwanted types
            activities = activities[-activities.ActivityTypeCode.isin(bad_activity_types)]
            activities = activities[-activities.ActivityMediaSubdivisionName.isin(bad_ActivityMediaSubdivisionNames)]    
            #all_ActivityMediaSubdivisionNames.append(activities.ActivityMediaSubdivisionName.unique())
            # restrict to my stations
            activities = activities[activities.MonitoringLocationIdentifier.isin(list(station_info_df.MonitoringLocationIdentifier.unique()))]
            activity_df = pd.concat([activity_df, activities])
                        
            #! why are there non-unique ActivityIdentifier?
            #! check other columns
            #date = activities.ActivityStartDate
            #activity_lst = activities.ActivityIdentifier
            #type_lst = activities.ActivityMediaSubdivisionName
            #activity_info_lst = list(map(lambda x, y:(x,y), activity_lst, type_lst))
            activity_info_df = activities[['ActivityIdentifier','ActivityMediaSubdivisionName']]            
   
            print('Loading activity for ' + str(m[0]) + " to " + str(m[1]) + " " + str(d))
        except:
            problems.append((m,d,'activity'))
            print('Could not load activity for ' + m[0] + " " + str(d))            
        try:
            result_url = url_prefix + "Result" + url_part_2 + str(m[0]) + url_post_date_01 + str(d) + url_post_date_02 + str(m[1]) + url_post_date_03 + str(d) + url_end + "&dataProfile=narrowResult"
            results = pd.read_csv(result_url, compression='zip', sep=',', low_memory=False)
            results = results.drop(bad_result_columns, axis=1)
            results = results[-results.CharacteristicName.isin(bad_CharacteristicNames)]
            # restrict to my activities
            results = results[results.ActivityIdentifier.isin(list(activity_info_df.ActivityIdentifier.unique()))]
            # add chem id
            results['chem_id'] = results.apply(lambda row: assign_chem_id(row), axis=1)
            results = results[results.chem_id.notnull()]
                # add bloq binary column
            results['bloq'] = results.apply(lambda row: is_bloq(row), axis=1)
            results = results[~results['ResultMeasure/MeasureUnitCode'].isin(bad_units)]            
            results['USGSPCode']= results['USGSPCode'].fillna(0.0).astype(int)           
            results = results[-results['USGSPCode'].isin([air_codes,
                                                          bed_sed_codes,
                                                          biota_codes,
                                                          solid_codes])]
            results['method_1'] = [1 if x in list(ss_codes) else 2 if x in list(dissolved_codes) else 3 if x in list(bulk_codes) else 0 for x in results.USGSPCode]
            results['method_2'] = [1 if x == 'Suspended' else 2 if x in ["Dissolved", 'Non-volatile', 'Volatile', 'Semivolatile'] else 3 if x in ['Total', 'Total Recoverable', 'Unfiltered'] else 0 for x in results.ResultSampleFractionText]
            results['method_3'] = [1 if x in [a[0] for a in ss_lst] else 2 if x in [b[0] for b in dissolved_lst] else 3 if x in [c[0] for c in bulk_lst] else 0 for x in results['ResultAnalyticalMethod/MethodIdentifier']]
            print('Loading results for ' + str(m[0]) + " to " + str(m[1]) + " " + str(d))
        except:
            problems.append((m,d,'result'))
            print('Could not load results for ' + m[0] + " " + str(d))

        try:
            loq_url = url_prefix + 'ResultDetectionQuantitationLimit' + url_part_2 + str(m[0]) + url_post_date_01 + str(d) + url_post_date_02 + str(m[1]) + url_post_date_03 + str(d) + url_end
            loq_df = pd.read_csv(loq_url, compression='zip', sep=',', low_memory=False)
            # select relevant columns
            loq_df = loq_df[['ResultIdentifier',
                             'CharacteristicName',
                             'DetectionQuantitationLimitTypeName',
                             'DetectionQuantitationLimitMeasure/MeasureValue',
                             'DetectionQuantitationLimitMeasure/MeasureUnitCode'
                             ]]
            #loq_df['loq_type'] = [1 if x in [a[0] for a in ss_lst] else 2 if x in [b[0] for b in dissolved_lst] else 3 if x in [c[0] for c in bulk_lst] else 0 for x in loq_df['DetectionQuantitationLimitTypeName']]
            loq_df['loq_type'] = [1 if x in technical_min else 2 if x in technical_qnt else 3 if x in reporting_min else 0 for x in loq_df['DetectionQuantitationLimitTypeName']]

            final_res = pd.merge(results,
                                 loq_df,
                                 on='ResultIdentifier',
                                 how='left')
            final_res['loq_val'] = final_res.apply(lambda row: calc_loq_val(row), axis=1)
            final_res['res_val'] = final_res.apply(lambda row: calc_res_val(row), axis=1)
            final_res['phase'] = final_res.apply(lambda row: get_method(row), axis=1)
            output_filename = 'result_' + m[0] + str(d) + '.csv'
            final_res.to_csv(output_filename)
            print('Loading loq for ' + str(m[0]) + " to " + str(m[1]) + " " + str(d))
        except:
            problems.append((m,d,'loq'))
            print('Could not load loq for ' + m[0] + " " + str(d))

#plot_df['loq_type'] = [1 if x in technical_min else 2 if x in technical_qnt else 3 if x in reporting_min else 0 for x in plot_df['DetectionQuantitationLimitTypeName']]

#Total and Total Recoverable same

# if DOC at the 'ActivityStartDate', 'MonitoringLocationIdentifier'
# if pH

#bad_ResultSampleFractionText = ['Bed Sediment',]

#bad_ResultStatusIdentifier  = []
   
#added opera val

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

plt.figure(figsize=(22,8))
ax = sns.boxplot(data = plot_method_w_pc[plot_method_w_pc.res_val > 0],
                hue = 'method', # different colors for different 'cls'
                x = 'DTXSID',
                y = 'ln_res_val')
ax2 = ax.twinx()
plot_method_w_pc.plot.scatter(x="DTXSID", y=np.log("WATER_SOLUBILITY_MOL/L_OPERA_PRED"), ax=ax2, legend=False, color="r")
ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
#sns.plt.title('Boxplot grouped by cls') # You can change the title here
plt.show()

import numpy as np
import statsmodels.api as sm # recommended import according to the docs
import matplotlib.pyplot as plt

sample = np.random.uniform(0, 1, 50)
ecdf = sm.distributions.ECDF(plot_df[plot_df.DTXSID == 'DTXSID4021395']['loq_val'])

x = np.linspace(min(sample), max(sample))
y = ecdf(x)
plt.step(x, y)
plt.show()

test_val = ['48HR_DAPHNIA_LC50_MOL/L_TEST_PRED',
            '96HR_FATHEAD_MINNOW_MOL/L_TEST_PRED',
            'TETRAHYMENA_PYRIFORMIS_IGC50_MOL/L_TEST_PRED']

plot_toxval['min_toxval'] = plot_toxval.loc[:, test_val].apply(pd.to_numeric, errors='coerce').min(axis=1)
plot_toxval['pnec'] = plot_toxval.min_toxval * plot_toxval.AVERAGE_MASS*1000
dtxsids_pcp[['chem_id','pnec']].drop_duplicates().to_csv('dtxsids_physchemprop_v2.csv')

plot_df = pd.DataFrame(columns = ['Activity Date','ResultIdentifier','CharacteristicName_x', 'DTXSID','bloq', 'loq_val','res_val','loq_type', 'method'])
for d in date_range:
    filename = 'filtered_result_'+ str(d) +'_final.csv'
    df = pd.read_csv(filename, low_memory=False)
    df['loq_type'] = [1 if x in technical_min else 2 if x in technical_qnt else 3 if x in reporting_min else 0 for x in df['DetectionQuantitationLimitTypeName']]
    df['method'] = df.apply(lambda row: get_method(row), axis=1)
    plot_df = pd.concat([plot_df, df[['ActivityStartDate','ResultIdentifier','CharacteristicName_x', 'DTXSID','bloq', 'loq_val','res_val','loq_type', 'method']]])
    

# =============================================================================
# for i in range(len(plot_loq.DTXSID.unique())):
#     this_chem = plot_loq.DTXSID.unique()[i]
#     this_set = plot_loq[plot_loq.DTXSID == this_chem]
#     try:
#         res_0_1 = ks_2samp(this_set[this_set.loq_type == 0].loq_val, this_set[this_set.loq_type == 1].loq_val)
#     except:
#         res_0_1 = np.nan
#     try:    
#         res_0_2 = ks_2samp(this_set[this_set.loq_type == 0].loq_val, this_set[this_set.loq_type == 2].loq_val)
#     except:
#         res_0_2 = np.nan
#     try:
#         res_0_3 = ks_2samp(this_set[this_set.loq_type == 0].loq_val, this_set[this_set.loq_type == 3].loq_val)
#     except:
#         res_0_3 = np.nan
#     try:
#         res_0_4 = ks_2samp(this_set[this_set.loq_type == 0].loq_val, this_set[this_set.loq_type == 4].loq_val)
#     except:
#         res_0_4 = np.nan
#     try:    
#         res_1_2 = ks_2samp(this_set[this_set.loq_type == 1].loq_val, this_set[this_set.loq_type == 2].loq_val)
#     except:
#         res_1_2 = np.nan
#     try:
#         res_1_3 = ks_2samp(this_set[this_set.loq_type == 1].loq_val, this_set[this_set.loq_type == 3].loq_val)
#     except:
#         res_1_3 = np.nan
#     try:    
#         res_1_4 = ks_2samp(this_set[this_set.loq_type == 1].loq_val, this_set[this_set.loq_type == 4].loq_val)
#     except:
#         res_1_4 = np.nan
#     try:
#         res_2_3 = ks_2samp(this_set[this_set.loq_type == 2].loq_val, this_set[this_set.loq_type == 3].loq_val)
#     except:
#         res_2_3 = np.nan
#     try:    
#         res_2_4 = ks_2samp(this_set[this_set.loq_type == 2].loq_val, this_set[this_set.loq_type == 4].loq_val)
#     except:
#         res_2_4 = np.nan
#     try:    
#         res_3_4 = ks_2samp(this_set[this_set.loq_type == 3].loq_val, this_set[this_set.loq_type == 4].loq_val)
#     except:
#         res_3_4 = np.nan
#     ks_lst.append((this_chem,res_0_1, res_0_2, res_0_3, res_0_4, res_1_2, res_1_3, res_1_4, res_2_3, res_2_4, res_3_4))
#     
# ntest_0_2 = len(this_set[this_set.loq_type == 2].loq_val)/len(this_set[this_set.loq_type == 0].loq_val)
# 
# =============================================================================

   

brief_result_columns = ['ResultIdentifier', #unique id
                        'ActivityStartDate', #for season
                        'MonitoringLocationIdentifier',                        
                        'loq_type', #for loq type
                        'loq_val',
                        'res_val', 
                        'bloq',
                        'method', #for phase
                        'chem_id']
small_df = pd.DataFrame(columns=brief_result_columns)

for d in dates:
    for m in months:
        filename = 'result_' + m[0] + str(d) + '_w_chem_id.csv'
        df = pd.read_csv(filename)
        small_df = pd.concat([small_df,df[brief_result_columns]],ignore_index='True')

def get_season(row):
    this_month = row['ActivityStartDate'][5:-3]
    
    if this_month in seasons[0]:
        return 3
    elif this_month in seasons[1]:
        return 0
    elif this_month in seasons[2]:
        return 1
    elif this_month in seasons[3]:
        return 2
    else: return np.na()

small_df['season'] = small_df.apply(lambda row: get_season(row), axis=1)
small_df['season'] = small_df['season'].astype(str)
small_df['phase'] = small_df['phase'].astype(str)
small_df['loq_type'] = small_df['loq_type'].fillna(0).astype(int)
small_df['loq_type'] = small_df['loq_type'].astype(int)

small_df['loq_type'] = small_df['loq_type'].astype(str)

#drop chem_id 0, both na
both_na = [167, #check these again
176,
178,
194,
195,
201,
206,
269,
362,
433,
442,
452,
462,
469,
470,
472,
538,
547,
562,
575,
632,
662,
824,
825,
831,
834,
839,
841,
856,
860,
865,
866,
874,
875,
880,
889,
891,
893,
894,
895,
897,
899,
900,
901,
902,
905,
907,
908,
915,
918,
920,
921,
922,
926,
927,
928,
929,
935,
936,
937,
938,
939,
945,
955,
958,
960,
966,
971,
972,
977,
978,
983,
984,
987,
991,
1001,
1005,
1019,
1028,
1029,
1043,
1046,
1049,
1051,
1062,
1065,
1067,
1069,
1071,
1533,
1555,
1557]

both_na = small_df[small_df['loq_val'].isna() & small_df['res_val'].isna()].groupby('chem_id')['chem_id'].count()
          
#join station info df LATER 
#! co-occurrence?

from scipy.stats import ks_2samp 
from itertools import combinations
from scipy.stats import wasserstein_distance

to_comb_lst = [0,1,2,3]
comb_lst = [",".join(map(str, comb)) for comb in combinations(to_comb_lst, 2)]
    

#http://sparky.rice.edu/astr360/kstest.pdf
#https://www.graphpad.com/guides/prism/7/statistics/interpreting_results_kolmogorov-smirnov_test.htm
ks_lst = []

#! A value is trying to be set on a copy of a slice from a DataFrame.
#Try using .loc[row_indexer,col_indexer] = value instead
def get_ks(df,val,category):
    for i in range(len(chem_id_df) - 1):
        this_chem_id = i + 1
        #this_chem_name = chem_id_df.CharacteristicName.iloc[i]
        this_set = df[df.chem_id == float(this_chem_id)]
        
        for c in comb_lst:
            n_1 = c[0]
            n_2 = c[2]
            try:            
                res = ks_2samp(this_set[this_set[category] == n_1][val].values, this_set[this_set[category] == n_2][val].values)
                ksstat= res[0]
                pval = res[1]
                ntest = len(this_set[this_set[category] == n_2][val])/len(this_set[this_set[category] == n_1][val]) 
                
                W = wasserstein_distance(this_set[this_set[category] == n_1][val].values, this_set[this_set[category] == n_2][val].values)
            except:
                res = np.nan
                ksstat = np.nan
                pval = np.nan
                ntest = np.nan
                W = np.nan
            ks_lst.append((this_chem_id, c, ksstat, pval, ntest, W))

get_ks(small_df, 'res_val', 'season')
ks_df = pd.DataFrame(ks_lst, columns=['this_chem_id', 'c', 'ksstat', 'pval', 'ntest', 'W'])
ks_df[ks_df['pval'] < 0.05].groupby('c')['this_chem_id'].count()/ks_df[ks_df['pval'] > 0.05].groupby('c')['this_chem_id'].count()

ax = small_df[small_df['res_val']>0]['season'].value_counts().plot(kind='bar',
                                    figsize=(14,8))

# success 'DTXSID2021319'

import matplotlib.pyplot as plt
import seaborn as sns
dtxsids_pcp = pd.read_csv('dtxsids_physchemprop.csv')

for_heatmap = dtxsids_pcp[['chem_id','AVERAGE_MASS','ATMOSPHERIC_HYDROXYLATION_RATE_(AOH)_CM3/MOLECULE*SEC_OPERA_PRED',
       'BIODEGRADATION_HALF_LIFE_DAYS_DAYS_OPERA_PRED',
       'BOILING_POINT_DEGC_OPERA_PRED', 'HENRYS_LAW_ATM-M3/MOLE_OPERA_PRED',
       'OCTANOL_AIR_PARTITION_COEFF_LOGKOA_OPERA_PRED',
       'OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED',
       'MELTING_POINT_DEGC_OPERA_PRED', 'VAPOR_PRESSURE_MMHG_OPERA_PRED',
       'WATER_SOLUBILITY_MOL/L_OPERA_PRED']].drop_duplicates()

for_heatmap.replace(to_replace='-',value=np.nan, inplace=True)
for_heatmap['Henrys Law ln(atm*m^3/mol)'] = np.log(for_heatmap['HENRYS_LAW_ATM-M3/MOLE_OPERA_PRED'].astype(float))
for_heatmap['Vapor pressure ln(mmHg)'] = np.log(for_heatmap['VAPOR_PRESSURE_MMHG_OPERA_PRED'].astype(float))
for_heatmap['log Octanol:air'] = np.log(for_heatmap['OCTANOL_AIR_PARTITION_COEFF_LOGKOA_OPERA_PRED'].astype(float))
for_heatmap['log Octanol:water'] = np.log(for_heatmap['OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED'].astype(float))
for_heatmap['Water solubility ln(mg/L)'] = np.log(for_heatmap['WATER_SOLUBILITY_MOL/L_OPERA_PRED'].astype(float) * for_heatmap['AVERAGE_MASS'].astype(float))

smaller_heatmap = for_heatmap[['Henrys Law ln(atm*m^3/mol)','log Octanol:air','log Octanol:water',
       'Water solubility ln(mg/L)']]


plt.figure(figsize=(20,20))
sns.heatmap(smaller_heatmap.astype(float).sort_values(by='WATER_SOLUBILITY_MOL/L_OPERA_PRED'), robust=True,cmap="YlGnBu")

sns.scatterplot(to_plot.chem_id, to_plot.ln_obs, alpha=0.3)

#plot 
plt.figure(figsize=(22,8))
ax = sns.boxplot(data = to_plot[to_plot.pnec.notnull()],
                #hue = 'method', # different colors for different 'cls'
                x = 'chem_id',
                y = 'ln_obs')
ax2 = ax.twinx()
plt.vlines(x=to_plot[to_plot.pnec.notnull()].chem_id, ymin=np.log(to_plot[to_plot.pnec.notnull()].mean_MLE), ymax=np.log(to_plot[to_plot.pnec.notnull()].mean_MLE+to_plot[to_plot.pnec.notnull()].sd_MLE))
to_plot.plot.scatter(x="chem_id", y="V3", ax=ax2, legend=False, color="r")
ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
#sns.plt.title('Boxplot grouped by cls') # You can change the title here
plt.show()

bar_heights = []
type_zero = []
type_one = []
type_two = []
type_three = []
for chem in pct_cen['chem_id'].unique():
    this_set = pct_cen[pct_cen.chem_id == chem]
    bar_height = len(this_set[this_set.bloq == 1])/len(this_set)
    type0 = len(this_set[(this_set.bloq == 1) & (this_set.loq_type == 0)])/len(this_set)
    bar_heights.append(bar_height)
    type_zero.append(type0)


N = len(pct_cen.chem_id.unique())
#menMeans = (20, 35, 30, 35, 27)
#womenMeans = (25, 32, 34, 20, 25)
#menStd = (2, 3, 4, 1, 2)
#womenStd = (3, 5, 2, 3, 3)
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence

p1 = plt.bar(ind, type_zero, width)
p2 = plt.bar(ind, type_one, width,
             bottom=type_zero)
p3 = plt.bar(ind, type_two, width,
             bottom=type_one)
p4 = plt.bar(ind, type_three, width,
             bottom=type_two)

plt.ylabel('Percent censored')
plt.title('Percent censored per chemical by limit type')
plt.xticks(ind, pct_cen['chem_id'].unique())
#plt.yticks(np.arange(0, 81, 10))
plt.legend((p1[0], p2[0], p3[0], p4[0]), ('Unknown', 'Tech min', 'Tech qnt','Rprt min'))

plt.show()

def set_chem_id(row):
    if row['all_chem_id'] > 0:
        return row['all_chem_id']
    elif row['bulk_chem_id'] > 0:
        return row['bulk_chem_id']
    elif row['dslv_chem_id'] > 0:
        return row['dslv_chem_id']
    else: return 0

loading_conc_df = pd.merge(loading_df[loading_df.DTXSID != "-"],abd_df,left_on='chem_id_y', right_on='all_chem_id',how='outer')

plt.figure(figsize=(22,8))
a = plt.scatter(x=plot_compare['chem_id'], y=plot_compare['ln_geom_mean_rmv'], c='black', size=5)
b = plt.scatter(x=plot_compare['chem_id'], y=plot_compare['ln_all_median_K-M'])
c = plt.scatter(x=plot_compare['chem_id'], y= plot_compare['ln_all_median_ROS'])
plt.scatter(x=plot_compare['chem_id'], y=plot_compare['ln_all_median_MLE'])
plt.scatter(x=plot_compare['chem_id'], y= plot_compare['ln_all_mean_MLE'])
plt.scatter(x=plot_compare['chem_id'], y= plot_compare['ln_all_mean_ROS'])
plt.scatter(x=plot_compare['chem_id'], y=plot_compare['ln_all_mean_K-M'])
plt.scatter(x=plot_compare['chem_id'], y= plot_compare['ln_bulk_median_K-M'], marker="s")
plt.scatter(x=plot_compare['chem_id'], y= plot_compare['ln_bulk_median_ROS'], marker="s")
plt.scatter(x=plot_compare['chem_id'], y=plot_compare['ln_bulk_median_MLE'], marker="s")
plt.scatter(x=plot_compare['chem_id'], y= plot_compare['ln_bulk_mean_MLE'], marker="s")
plt.scatter(x=plot_compare['chem_id'], y= plot_compare['ln_bulk_mean_ROS'], marker="s")
plt.scatter(x=plot_compare['chem_id'], y=plot_compare['ln_bulk_mean_K-M'], marker="s")
plt.scatter(x=plot_compare['chem_id'], y= plot_compare['ln_dslv_median_K-M'], marker="^")
plt.scatter(x=plot_compare['chem_id'], y= plot_compare['ln_dslv_median_ROS'], marker="^")
plt.scatter(x=plot_compare['chem_id'], y=plot_compare['ln_dslv_median_MLE'], marker="^")
plt.scatter(x=plot_compare['chem_id'], y= plot_compare['ln_dslv_mean_MLE'], marker="^")
plt.scatter(x=plot_compare['chem_id'], y= plot_compare['ln_dslv_mean_ROS'], marker="^")
plt.scatter(x=plot_compare['chem_id'], y=plot_compare['ln_dslv_mean_K-M'], marker="^")
plt.legend()


plt.scatter(plot_load['dslv_mean_K-M'], plot_load['dslv_median_K-M'], color='b', marker="^")
plt.scatter(plot_load['bulk_mean_K-M'], plot_load['bulk_median_K-M'], color='b', marker="s")
plt.yscale('log')
plt.xscale('log')

dtxsids_pcp['OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED'] = pd.to_numeric(dtxsids_pcp['OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED'], errors='coerce')