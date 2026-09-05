#!/bin/bash
# sequential queue after F3 finishes
cd /home/user/H16P/audit/fable_engine
while pgrep -f "sweep_shi.py data/F3_lam0_L4" > /dev/null; do sleep 30; done
python3 evolve.py kkl data/F4_kkl_evolve.jsonl 40 96 0.05 --dims=3,4,9,10,11,7,8 > data/F4_kkl_evolve.log 2>&1
python3 evolve.py yz data/F8_yz_evolve.jsonl 40 96 0.05 > data/F8_yz_evolve.log 2>&1
python3 sweep_family.py q4pert data/F7_q4pert.jsonl 20000 5 --eps=1e-3 > data/F7_q4pert.log 2>&1
python3 sweep_family.py mv data/F9_mv.jsonl 40000 7 --L=1.0 > data/F9_mv.log 2>&1
python3 sweep_family.py kklx data/F5_kklx.jsonl 60000 9 --L=1.5 > data/F5_kklx.log 2>&1
python3 sweep_shi.py data/F5_shi_full_L4.jsonl 150000 13 --L=4 > data/F5_shi_full_L4.log 2>&1
python3 evolve.py shi data/F8_shi_evolve.jsonl 30 96 0.05 > data/F8_shi_evolve.log 2>&1
python3 sweep_family.py q3rpert data/F6_q3rpert.jsonl 20000 17 --eps=1e-3 > data/F6_q3rpert.log 2>&1
python3 sweep_family.py mvpert data/F9_mvpert.jsonl 20000 19 --eps=1e-3 > data/F9_mvpert.log 2>&1
python3 sweep_shi.py data/F3b_lam0_store1.jsonl 60000 23 --lam0 --L=4 --store1 > data/F3b.log 2>&1
python3 evolve.py file:data/F3b_lam0_store1.jsonl data/F3c_lam0_evolve.jsonl 40 96 0.1 --dims=3,4,5,9,10 > data/F3c.log 2>&1
echo QUEUE DONE > data/QUEUE_DONE
