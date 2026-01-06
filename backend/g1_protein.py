import requests
from collections import Counter
import networkx as nx
import matplotlib.pyplot as plt 

def protein_summary(protein_id:str):
    try:
        uniprot_id = protein_id

        # UniProt REST API URL (JSON format)
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"

        response = requests.get(url)
        data = response.json()
        
        protein_name = data["proteinDescription"]["recommendedName"]["fullName"]["value"]
        sequence = data["sequence"]["value"]
        if not protein_name or not sequence:
            return None
        

        counts = Counter(sequence)

        most_common = counts.most_common(1)[0]
        least_common = counts.most_common()[-1] 

        # domains
        url_domain = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{protein_id}/"
        response = requests.get(url_domain).json()
        domain = []
        for item in response["results"]:
            meta = item["metadata"]
            domain.append((meta["accession"],meta.get("name")))
        domain_text = "\n".join([f"{a} - {n}" for (a, n) in domain])
        
        result = (
            f"Protein Code: {protein_id}\n"
            f"Predicted Name: {protein_name}\n"
            f"Sequence Length: {len(sequence)} amino acids\n"
            f"Amino Acid Sequence: {sequence})\n"
            f"Most Frequent Amino Acid: {most_common[0]} ({most_common[1]} times)\n"
            f"Least Frequent Amino Acid: {least_common[0]} ({least_common[1]} times)\n"
            f"Domains:\n{domain_text}"
        )
        return result
    except Exception:
        return None

## ----------------- Function 2 ---------------------------------

def get_alphafold_pdb(protein_id):
    """
    Returns PDB string from AlphaFold
    """
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{protein_id}"
    response = requests.get(url)
    if response.status_code != 200:
        return None
    model = response.json()[0]
    pdb_url = model.get("pdbUrl")
    if not pdb_url:
        return None
    pdb_data = requests.get(pdb_url).text

    alpha = (
        f"=== ALPHAFOLD STRUCTURE INFORMATION ===\n"
        f"UniProt ID: {protein_id}\n"
        f"Model Confidence (pLDDT): {model.get('plddt', 'N/A')}\n"
        f"\nDownload Links:\n"
        f"PDB File : {model.get('pdbUrl', 'Not available')}\n"
        f"CIF File : {model.get('cifUrl', 'Not available')}\n"
        f"JSON Info: {url}\n"
    )
    return pdb_data, alpha



## ------------- Function 3 ----------------------------------
def ppi_network(protein_id):
        
    uniprot_id = protein_id

    url = f"https://string-db.org/api/json/network?identifiers={uniprot_id}&species=9606"

    response = requests.get(url)
    int_list = []

    if response.status_code != 200:
        return KeyError
    else:
        data = response.json()
   
        if len(data) == 0:
            return ("No interactions found for this protein.")
        else:
            #print("=== PROTEIN-PROTEIN INTERACTIONS (STRING DB) ===\n")
            for interaction in data:
                partner = interaction.get("preferredName_B", "N/A")
                score = interaction.get("score", "N/A")
                int_list.append((partner, score))

    G = nx.Graph()
    G.add_node(protein_id)

    for partner, score in int_list:
        G.add_node(partner)
        G.add_edge(protein_id, partner, weight=score)

    fig, ax = plt.subplots(figsize=(6, 5))
    fig.set_facecolor("#e7f2ff") 
    ax.set_facecolor("#e7f2ff")

    pos = nx.spring_layout(G, k=0.6, seed=42)

    nx.draw_networkx_nodes(G, pos, node_size=900, node_color="#5aa3e8", ax=ax)

    weights = [G[u][v]['weight'] * 3 for u, v in G.edges()] #thickness
    nx.draw_networkx_edges(G, pos, width=weights, edge_color="#2c6faa", ax=ax)

    nx.draw_networkx_labels(G, pos, font_size=9, font_color="white", ax=ax)

    ax.set_title(f"Protein Interaction Network ({protein_id})", fontsize=12)
    ax.axis("off")

    explain = "This network shows the predicted protein–protein interactions for your protein of interest. Each node represents a protein, and each edge represents an interaction. The thickness/strength of edges indicates the confidence level of the interaction. This visualization helps identify functional partners and potential pathways involving the protein."

    return fig, explain