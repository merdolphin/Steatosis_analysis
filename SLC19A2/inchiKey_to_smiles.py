import requests

# List of InChIKeys
#inchikeys = [
#    "NFGIENSPALNOON-UHFFFAOYSA-N",
#    "SEOVTRFCIGRIMH-UHFFFAOYSA-N",
#    "VVHOUVWJCQOYGG-REOHCLBHSA-N",
    # Add more InChIKeys here
#]

# Base URL for PubChem REST API
base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# Function to retrieve SMILES string for a given InChIKey
def get_smiles_for_inchikey(inchikey):
    # URL for InChIKey to SMILES conversion
    url = f"{base_url}/compound/inchikey/{inchikey}/property/CanonicalSMILES/JSON"
    
    try:
        # Make GET request to PubChem API
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for 4XX and 5XX status codes
        data = response.json()
        
        # Extract SMILES string from the response
        smiles = data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
        return smiles
    except Exception as e:
        print(f"Failed to retrieve SMILES for {inchikey}: {e}")
        return None


with open("inchiKeys", "r") as file:
    inchikeys = [line.strip() for line in file]


# Loop through each InChIKey and retrieve SMILES
for inchikey in inchikeys:
    smiles = get_smiles_for_inchikey(inchikey)
    if smiles:
        print(f"InChIKey: {inchikey}, SMILES: {smiles}")

