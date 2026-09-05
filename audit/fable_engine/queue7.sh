#!/bin/bash
cd /home/user/H16P/audit/fable_engine
python3 recount_hits.py data/F15_kklstar_sweep.jsonl data/F15_mvneutral_sweep.jsonl data/F18c_scaled.jsonl data/F21_shi_compact.jsonl > data/RECOUNT_fixed_counter.log 2>&1
python3 f19_q4_alien.py > data/F19_q4_alien_fixed.log 2>&1
echo DONE > data/QUEUE7_DONE
