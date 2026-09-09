"""Endpoint continuation from previously fitted anchor shifts +/-0.20."""
import json
from q4_anchor_placement import main,P
records=json.loads((P/'q4_anchor_placement.json').read_text())['records']
warm={(q['r'],1 if q['shift']>0 else -1):q['normalized_controls'] for q in records if abs(q['shift'])==.2}
main(amounts=(.225,.24,.2475,.249,.2499),outname='q4_anchor_limits.json',warm=warm)
