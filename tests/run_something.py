import pathlib

import arrow
import pandas as pd
from tqdm import tqdm
import selfies as sf
from selfies.utils.smiles_utils import tokenize_smiles
from selfies_rust.selfies_rust import encoder

TEST_DIR = pathlib.Path(__file__).parent
TEST_SET_PATH = TEST_DIR / "test_sets"

df = pd.read_csv(f"{TEST_SET_PATH}/zinc.csv")

start_time = arrow.now().timestamp()
for smi in df["smiles"]:
    encoder(smi, None, None)
print(arrow.now().timestamp() - start_time)

start_time = arrow.now().timestamp()
for smi in df["smiles"]:
    tokens = []
    for token in tokenize_smiles(smi):
        tokens.append(token)
print(arrow.now().timestamp() - start_time)