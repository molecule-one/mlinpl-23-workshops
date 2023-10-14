# NOTE: you might need to first restart database and remove solution_ folders.

#for protein in GSK3β DRD2_server JNK3; do
#  PYTHONPATH=`pwd` PORT=8000 python solutions/run.py -b 1000 -w ml -t $protein -u baseline-ml
#  PYTHONPATH=`pwd` PORT=8000 python solutions/run.py -b 1000 -w mutate -t $protein -u baseline-mutate
#  PYTHONPATH=`pwd` PORT=8000 python solutions/run.py -b 1000 -w mutate_low_exp -t $protein -u baseline-mutate-low-exp
#  PYTHONPATH=`pwd` PORT=8000 python solutions/run.py -b 1000 -w random -t $protein -u baseline-random
#done

for protein in DRD2_server GSK3β JNK3; do
  PYTHONPATH=`pwd` PORT=8000 python solutions/run.py -b 5000 -w mutate -t $protein -u baseline-mutate -s 25
  PYTHONPATH=`pwd` PORT=8000 python solutions/run.py -b 5000 -w mutate_low_exp -t $protein -u baseline-mutate-low-exp -s 25
  PYTHONPATH=`pwd` PORT=8000 python solutions/run.py -b 5000 -w ml -t $protein -u baseline-ml -s 25
  PYTHONPATH=`pwd` PORT=8000 python solutions/run.py -b 5000 -w random -t $protein -u baseline-random -s 25
done

#for protein in GSK3β DRD2_server JNK3; do
#  PYTHONPATH=`pwd` PORT=8000 python solutions/run.py -b 5000 -w ml -t $protein -u baseline-ml -s 50
#  PYTHONPATH=`pwd` PORT=8000 python solutions/run.py -b 5000 -w mutate -t $protein -u baseline-mutate -s 50
#  PYTHONPATH=`pwd` PORT=8000 python solutions/run.py -b 5000 -w mutate_low_exp -t $protein -u baseline-mutate-low-exp -s 50
#  PYTHONPATH=`pwd` PORT=8000 python solutions/run.py -b 5000 -w random -t $protein -u baseline-random -s 50
#done