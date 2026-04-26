import urllib.request
import json
import time

def fetch_smiles_by_name(name: str):
    """
    Fetches the SMILES string for a given chemical name from PubChem.
    """
    print(f"[*] Fetching SMILES for '{name}' from PubChem...")
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(name)}/property/CanonicalSMILES,IsomericSMILES/JSON"
        with urllib.request.urlopen(url, timeout=10) as response:
            data = json.loads(response.read().decode())
            props = data['PropertyTable']['Properties'][0]
            # Try various keys PubChem might use
            for key in ['CanonicalSMILES', 'IsomericSMILES', 'SMILES', 'ConnectivitySMILES']:
                if key in props:
                    smiles = props[key]
                    print(f"[*] Found SMILES for {name}: {smiles}")
                    return smiles
            print(f"[!] No SMILES found in PubChem response for '{name}'")
            return None
    except Exception as e:
        print(f"[!] Error fetching SMILES for '{name}': {e}")
        return None
