import os
from amis import reconstruct_multi_models
from Bio import SeqIO

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(
        description='Recommend substitutions to a wildtype sequence'
    )
    parser.add_argument(
        '--model-names',
        type=str,
        default=[ 'esm1b', 'esm1v1', 'esm1v2', 'esm1v3', 'esm1v4', 'esm1v5', ],
        nargs='+',
        help='Type of language model (e.g., esm1b, esm1v1)'
    )
    parser.add_argument(
        '--alpha',
        type=float,
        default=None,
        help='alpha stringency parameter'
    )
    parser.add_argument(
        '--cuda',
        type=str,
        default='cuda',
        help='cuda device to use'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='output file'
    )
    parser.add_argument(
        '--fasta',
        type=str,
        default=None,
        help='input fasta file'
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()

    if ":" in args.cuda:
        os.environ["CUDA_VISIBLE_DEVICES"] = args.cuda.split(':')[-1]

    fasta_file = args.fasta
    output_file = args.output

    with open(output_file, 'w') as f:
        print(f'pid,mutation,number_of_models', file=f)
        for record in SeqIO.parse(open(fasta_file), "fasta"):
            pid = record.id
            mutations_models = reconstruct_multi_models(
                record.seq,
                args.model_names,
                alpha=args.alpha,
            )
            for k, v in sorted(mutations_models.items(), key=lambda item: -item[1]):
                mut_str = f'{k[1]}{k[0] + 1}{k[2]}'
                print(f'{pid},{mut_str},{v}', file=f)
