#!/bin/bash
cd /home/user/H16P/audit/fable_engine
python3 sweep_log.py mvneutral data/F15_mvneutral_sweep.jsonl 20000 41 --store=3 > data/F15_mvneutral_sweep.log 2>&1
python3 sweep_log.py kklstar data/F15_kklstar_sweep.jsonl 10000 43 --store=3 > data/F15_kklstar_sweep.log 2>&1
echo DONE > data/F15_QUEUE_DONE
