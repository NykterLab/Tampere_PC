import argparse as ap
import datetime
import os

parser = ap.ArgumentParser(description='Train a BPNET model')
parser.add_argument('--dataspec', dest='dataspec', help='dataspec file name')
parser.add_argument('--model_dir', dest='model_dir', help='model directory to be created')
parser.add_argument('--configgin', dest='configgin', help='gin configuration file name')
parser.add_argument('--workers', dest='workers', help='number of workers')
args = parser.parse_args()
dataspec = args.dataspec
model_dir = args.model_dir
configgin = args.configgin
workers = args.workers

run_id = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

command = f"bpnet train {dataspec} {model_dir}\
                --premade=bpnet9 \
                --config=\"{configgin}\" \
                --run-id '{run_id}'\
                --num-workers {workers} \
                --in-memory \
                --vmtouch"

os.system(command)

os.symlink(f'{run_id}/seq_model.pkl', f'{model_dir}/seq_model.pkl')
