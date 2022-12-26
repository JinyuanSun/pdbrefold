import pandas as pd
import requests
import shutil
import gzip
import py3Dmol
from stmol import showmol
import streamlit as st
import numpy as np
from tmtools import tm_align


def download_file(url):
    local_filename = url.split('/')[-1]
    with requests.get(url, stream=True) as r:
        with open(local_filename, 'wb') as f:
            shutil.copyfileobj(r.raw, f)

    return local_filename


def fetch_pdb(uni_id):
    response = requests.get(
        f'https://alphafold.ebi.ac.uk/files/AF-{uni_id}-F1-model_v4.pdb')
    return response.content.decode('utf-8')




def update_mapping():
    # local_csv = download_file('ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz')
    with gzip.open('pdb_chain_uniprot.csv.gz', 'rb') as f:
        df = pd.read_csv(f, skiprows=1, sep=',')
    return df


def render_mol(pdb: str, start=0, stop=5):
    # pdb = open(pdb_path, "r").read()
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb, "pdb")
    pdbview.setStyle({"cartoon": {"color": "#d5f8fa"}})
    pdbview.setBackgroundColor("#0E1117")

    i = 0
    j = 0
    for line in pdb.split("\n"):
        split = line.split()
        color = "#d5f8fa"

        if len(split) == 0 or split[0] != "ATOM":
            continue
        if "CA" in split:
            j += 1
        if j >= stop:
            color = "grey"
        if j <= start:
            color = "grey"

        pdbview.setStyle({"model": -1, "serial": i},
                         {"cartoon": {"color": color}})
        i += 1

    pdbview.zoomTo()
    pdbview.zoom(2, 600)
    showmol(pdbview, height=500, width=600)


def get_pdb_subdf(pdb, chain):
    pdb_df = df[df['PDB'] == pdb.lower()]
    odf = pdb_df[pdb_df['CHAIN'] == chain]
    return odf


def get_afpdb(uniprot_id):
    afpdb_name = download_file(
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb")
    return afpdb_name


def get_pdb_chain(pdb_id, chain_id):
    response = requests.get(f'https://files.rcsb.org/download/{pdb_id}.pdb')
    res = ''
    for line in response.content.decode('utf-8').split("\n"):
        if line.startswith("ATOM"):
            if line[21] == chain_id:
                res += line+"\n"
    return res


def get_ca_coords_from_string(pdb_string):
    coords = []
    for line in pdb_string.split("\n"):
        if line.startswith("ATOM") and line[12:16].replace(" ", "") == 'CA':
            x = float(line[30: 38])
            y = float(line[38: 46])
            z = float(line[46: 54])
            coords.append([x, y, z])
    return np.array(coords)


if __name__ == '__main__':
    from tqdm import tqdm
    st.set_page_config("PDB-REFOLD", page_icon="🔖", layout="wide")
    df = update_mapping()

    st.sidebar.header('PDB-REFOLD')
    st.sidebar.write("Refold PDB structures with AlphaFold")
    x = st.sidebar.text_input("PDB ID (e.g.: *5xjh*)", value="5xjh")
    y = st.sidebar.text_input("Chain ID (e.g.: *A*)", value="A")
    subdf = get_pdb_subdf(x, y)
    # print(subdf)
    uniprot_id = subdf['SP_PRIMARY'].values[0]
    start = subdf['SP_BEG'].values[0]
    stop = subdf['SP_END'].values[0]
    # print(uniprot_id, start, stop)

    af_string = fetch_pdb(uniprot_id)
    pdb_string = get_pdb_chain(x, y)
    pdb_coords = get_ca_coords_from_string(pdb_string)
    af_coords = get_ca_coords_from_string(af_string)
    res = tm_align(pdb_coords, af_coords, 'A' *
                   len(pdb_coords), 'A'*len(af_coords))
    print(res.tm_norm_chain2)

    col1, col2 = st.columns(2)
    with col1:
        render_mol(af_string, start, stop)
        st.write('AlpahFold Predicted')
    with col2:
        render_mol(pdb_string, -1, 9999)
        st.write('Crystal structure')
    st.write('### TM-Score normlized by crystal structure: ' +
             "{:.4f}".format(res.tm_norm_chain2))
    st.write(
        """## Citations\n
1. Jumper, J et al. Highly accurate protein structure prediction with AlphaFold. Nature (2021).\n
2. Varadi, M et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. Nucleic Acids Research (2021).\n
3. Rego, N et al. 3Dmol.js: molecular visualization with WebGL. Bioinformatics (2014).
    """
    )
    st.sidebar.download_button(
        'Use enter to update view and click this button to download', file_name=f"AF_{x}_{y}.pdb", data=pdb_string)
