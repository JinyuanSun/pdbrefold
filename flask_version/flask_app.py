import os
import pandas as pd
import requests
import shutil
import gzip
import numpy as np
from flask import Flask, render_template, request
from flask import redirect

# from flask_nglview import NGLView
from tmtools import tm_align

app = Flask(__name__)
# app.register_blueprint(NGLView)

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
    with gzip.open('pdb_chain_uniprot.csv.gz', 'rb') as f:
        df = pd.read_csv(f, skiprows=1, sep=',')
    return df

def get_pdb_subdf(pdb, chain):
    pdb_df = df[df['PDB'] == pdb.lower()]
    odf = pdb_df[pdb_df['CHAIN'] == chain]
    return odf

def get_pdb_chain(pdb_id, chain_id):
    response = requests.get(f'https://files.rcsb.org/download/{pdb_id}.pdb')
    res = ''
    for line in response.content.decode('utf-8').split("\n"):
        if line.startswith("ATOM"):
            if line[21] == chain_id:
                res += line+"\n"
    return res

def get_ca_coords_and_atom_lines_from_string(pdb_string):
    ca_coords = []
    ca_lines = []
    atom_lines = []
    for line in pdb_string.split("\n"):
        if line.startswith("ATOM"):
            atom_lines.append(line)
            if line[12:16].replace(" ", "") == 'CA':
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ca_coords.append([x, y, z])
                ca_lines.append(line)
    return np.array(ca_coords), ca_lines, atom_lines

def apply_tm_transform(coord, t, u):
    return np.dot(u, coord) + t

df = update_mapping()

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        pdb_id = request.form.get("pdb_id")
        chain_id = request.form.get("chain_id")
        subdf = get_pdb_subdf(pdb_id, chain_id)
        uniprot_id = subdf['SP_PRIMARY'].values[0]
        start = subdf['SP_BEG'].values[0]
        stop = subdf['SP_END'].values[0]

        af_string = fetch_pdb(uniprot_id)
        pdb_string = get_pdb_chain(pdb_id, chain_id)
        pdb_coords, pdb_ca_lines, pdb_atom_lines = get_ca_coords_and_atom_lines_from_string(pdb_string)
        af_coords, af_ca_lines, _ = get_ca_coords_and_atom_lines_from_string(af_string)
        res = tm_align(pdb_coords, af_coords, 'A' * len(pdb_coords), 'A' * len(af_coords))

        updated_pdb_lines = []
        for line in pdb_atom_lines:
                coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                updated_coord = apply_tm_transform(coord, res.t, res.u)
                updated_line = line[:30] + "{:8.3f}{:8.3f}{:8.3f}".format(*updated_coord) + line[54:]
                updated_pdb_lines.append(updated_line)

        updated_pdb_string = "\n".join(updated_pdb_lines)

        return render_template("results.html", af_string=af_string, pdb_string=updated_pdb_string, tm_score=res.tm_norm_chain2, uniprot_id=uniprot_id)

    return render_template("index.html")

# @app.route('/download_pdb/<pdb_id>/<chain_id>')
# def download_pdb(pdb_id, chain_id):
#     pdb_string = get_pdb_chain(pdb_id, chain_id)
#     response = app.response_class(
#     response=pdb_string,
#     status=200,
#     mimetype='application/octet-stream',
#     headers={"Content-Disposition": f"attachment;filename=AF_{pdb_id}_{chain_id}.pdb"}
#     )
#     return response

@app.route('/download_pdb/<uniprot_id>')
def download_pdb(uniprot_id):
    alphafold_db_url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'
    return redirect(alphafold_db_url)

if __name__ == "__main__":
    app.run(debug=True)
