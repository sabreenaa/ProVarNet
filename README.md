# ProVarNet

Protein Structure, Variant & Network Analysis Tool
ProVarNet is a PyQt5-based desktop application for exploring protein structure, sequence variants, disease associations, and protein–protein interaction networks using a single UniProt ID as input. It integrates multiple analysis modules into an interactive graphical interface designed for biological data exploration and interpretation.

✨ Features

Protein Summary Retrieval

  - Fetches and displays protein metadata from UniProt.

  - Clean and readable summary view.

AlphaFold Structure Viewer

  - Interactive 3D visualisation of protein structures using py3Dmol.

  - Colored structure display.

Protein–Protein Interaction Network

  - Visualises interaction networks as plots.

Variant Analysis

  - Graphical analysis of protein variants.
  
  - Summary and explanation panels for interpretation.

Disease-associated Variant Analysis

  - Summarise disease associations and visualise variant distributions.

  - Displays disease-linked variants in a sortable, filterable table.

  - Column-wise filtering.

## 🧬 Input -
The application accepts a UniProt protein ID (e.g. P04637, P38398, Q9Y6K9) as input.

# UI Overview

### Welcome Page

  - Input UniProt ID.

  - Fetch protein summary.

### Menu Page

Navigate to:

  - AlphaFold Structure

  - Protein–Protein Interactions

  - Variant Analysis

  - Disease Association

### Dialog Windows

Dedicated windows for each analysis:

  - Structure viewer

  - Interaction network

  - Variant distribution plots

  - Disease-associated variants table


⚠️ Notes

Requires an internet connection for API-based data retrieval.

Some proteins may not have variant or disease data available.

Large proteins may take longer to visualise.



## INSTALLATION & RUNNING THE APP

OPTION 1 (RECOMMENDED): DOWNLOAD ZIP

On the GitHub repository page.

Click "Code" and then "Download ZIP".

Extract (unzip) the folder.

Open a terminal or command prompt inside the extracted folder (where app.py is located).

### Run the app

```bash
pip install -r requirements.txt
python app.py
```

OPTION 2 (ADVANCED USERS): CLONE THE REPOSITORY
```bash
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
pip install -r requirements.txt
python app.py
```
OPTIONAL (RECOMMENDED): USE A VIRTUAL ENVIRONMENT

Create a virtual environment:

```bash
python -m venv venv
```
Activate it:

On Linux or macOS:
```bash
source venv/bin/activate
```
On Windows:

```bash
venv\Scripts\activate
```
Then install dependencies and run:
```bash
pip install -r requirements.txt
python app.py
```
TROUBLESHOOTING

Make sure you run:

python app.py

from the folder that contains app.py.

If PyQtWebEngine fails to install:

```bash
pip install pyqtwebengine
```

