import numpy as np
import sys
import json
from scivae import VAE


def main(work_matrix, config):
    vae = VAE(work_matrix, work_matrix, ["None"] * len(work_matrix), config, f'vae_rcm')
    vae.encode('default', epochs=20, batch_size=50, train_percent=75.0)
    return vae

if __name__ == '__main__':
    work_matrix = np.load(sys.argv[2])
    with open(sys.argv[1]) as js:
        config = json.load(js)
    prefix = sys.argv[3]

    vae = main(work_matrix, config)
    
    vae.save(weight_file_path=f'{prefix}_weights.h5', optimizer_file_path=f'{prefix}_optimiser.json',
             config_json=f'{prefix}_config.json')

    encoding = vae.get_encoded_data()
    np.save(f'{prefix}.embedding.npy', encoding)