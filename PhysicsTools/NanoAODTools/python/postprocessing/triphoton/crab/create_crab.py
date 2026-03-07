"""
Generate per-sample CRAB configs from sample JSON files.

Usage:
    python create_crab.py 2018
    python create_crab.py 2017
    python create_crab.py 2016
    python create_crab.py 2016apv
    python create_crab.py 2022
    python create_crab.py 2022EE
"""
import os
import sys
import json
from shutil import copyfile

year = sys.argv[1]

LUMI_JSONS = {
    "2016apv": "Cert_271036-284044_13TeV_Legacy2016_Collisions16_preVPF_JSON.txt",
    "2016": "Cert_271036-284044_13TeV_Legacy2016_Collisions16_postVPF_JSON.txt",
    "2017": "Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
    "2018": "Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
    "2022": "Cert_Collisions2022_355100_362760_Golden.json",
    "2022EE": "Cert_Collisions2022_355100_362760_Golden.json",
}

if year not in LUMI_JSONS:
    print("ERROR: Unknown year '%s'. Valid: %s" % (year, list(LUMI_JSONS.keys())))
    sys.exit(1)

workdir = "config_crab_%s/" % year
if not os.path.isdir(workdir):
    os.makedirs(workdir)

datajson = LUMI_JSONS[year]
samplejson = os.path.join("samples", "samples%s.json" % year)
scriptpath = "%s_script" % year

with open(samplejson, 'r') as fin:
    lines = json.loads(fin.read())

for key, value in lines.items():
    # Data samples: [OutputName, DAS_path, crab_script.sh]
    if len(value) == 3:
        copyfile('data_cfg.py', workdir + key + '_cfg.py')
        das_escaped = value[1].replace('/', 'sss')
        os.system(r'sed -i "6s/dummy/%s/g" %s' % (value[0], workdir + key + '_cfg.py'))
        if any(('_' + era) in value[0] for era in 'ABCDEFGH'):
            os.system(r'sed -i "13s/dummy/%s/g" %s' % (scriptpath + 'sss' + value[2], workdir + key + '_cfg.py'))
        os.system(r'sed -i "13s/sss/\//g" %s' % (workdir + key + '_cfg.py'))
        os.system(r'sed -i "15s/dummy/%s/g" %s' % (datajson, workdir + key + '_cfg.py'))
        os.system(r'sed -i "19s/dummy/%s/g" %s' % (das_escaped, workdir + key + '_cfg.py'))
        os.system(r'sed -i "19s/sss/\//g" %s' % (workdir + key + '_cfg.py'))
        os.system(r'sed -i "23s/dummy/%s/g" %s' % (datajson, workdir + key + '_cfg.py'))
        os.system(r'sed -i "26s/dummy/%s/g" %s' % (value[0], workdir + key + '_cfg.py'))

    # MC samples: [OutputName, DAS_path]
    if len(value) == 2:
        print(value)
        copyfile('mc_cfg.py', workdir + key + '_cfg.py')
        das_escaped = value[1].replace('/', 'sss')
        os.system(r'sed -i "6s/dummy/%s/g" %s' % (value[0], workdir + key + '_cfg.py'))
        os.system(r'sed -i "13s/dummy/%s/g" %s' % (scriptpath + 'sss' + 'crab_script.sh', workdir + key + '_cfg.py'))
        os.system(r'sed -i "13s/sss/\//g" %s' % (workdir + key + '_cfg.py'))
        os.system(r'sed -i "19s/dummy/%s/g" %s' % (das_escaped, workdir + key + '_cfg.py'))
        os.system(r'sed -i "19s/sss/\//g" %s' % (workdir + key + '_cfg.py'))
        os.system(r'sed -i "26s/dummy/%s/g" %s' % (value[0], workdir + key + '_cfg.py'))
