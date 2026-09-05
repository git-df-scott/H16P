#!/bin/bash
cd /home/user/H16P/audit/fable_engine
python3 evolve.py file:data/kklx_seeds.jsonl data/F14_kklx_seeds_evolve3.jsonl 60 96 0.03 > data/F14_kklx_seeds_evolve3.log 2>&1
python3 evolve.py yz data/F8_yz_evolve.jsonl 40 96 0.05 > data/F8_yz_evolve.log 2>&1
python3 evolve.py shi data/F8_shi_evolve.jsonl 30 96 0.05 > data/F8_shi_evolve.log 2>&1
python3 sweep_family.py q3rpert data/F6_q3rpert.jsonl 20000 17 --eps=1e-3 > data/F6_q3rpert.log 2>&1
python3 sweep_family.py mvpert data/F9_mvpert.jsonl 20000 19 --eps=1e-3 > data/F9_mvpert.log 2>&1
python3 sweep_family.py q4pert data/F7_q4pert.jsonl 20000 5 --eps=1e-3 > data/F7_q4pert.log 2>&1
python3 sweep_family.py kklx data/F5_kklx.jsonl 60000 9 --L=3.0 > data/F5_kklx.log 2>&1
python3 sweep_shi.py data/F3b_lam0_store1.jsonl 60000 23 --lam0 --L=4 --store1 > data/F3b.log 2>&1
python3 evolve.py file:data/F3b_lam0_store1.jsonl data/F3c_lam0_evolve.jsonl 40 96 0.1 --dims=3,4,5,9,10 > data/F3c.log 2>&1
echo QUEUE DONE > data/QUEUE_DONE
