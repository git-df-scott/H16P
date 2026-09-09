from pathlib import Path
import json
from q4_finite_continuation import main
if __name__=='__main__':
    source=json.loads((Path(__file__).resolve().parent/'q4_finite_continuation.json').read_text())
    warm={row['r']:row['normalized_controls'] for row in source['records'] if row['epsilon']==.4}
    main(epsvalues=(.7,1.,1.5,2.,3.,4.,6.,8.),outname='q4_finite_extended.json',warm=warm)
