#!/bin/bash
cd /home/user/H16P/audit/fable_engine
while pgrep -f "^python3 sweep_log.py kklstar" > /dev/null; do sleep 30; done
python3 evolve_log.py data/F15_kklstar_sweep.jsonl data/F20_kklstar_evolve.jsonl 40 64 0.03 > data/F20_kklstar_evolve.log 2>&1
echo DONE > data/F20_DONE
