import argparse
import multiprocessing

from local_conservation_analysis_pipeline import conservation_pipeline

N_CORES = round(multiprocessing.cpu_count()/2)


if __name__ == '__main__':
    config_file = './params.yaml'
    parser = argparse.ArgumentParser(
        description="Run conservation analysis pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        metavar="<file>",
        required=True,
        help="""path to config file in yaml format""",
    )
    parser.add_argument(
        "-n",
        "--n_cores",
        type=int,
        metavar="<int>",
        default=N_CORES,
        help=f"""number of cores to use""",
    )
    args = parser.parse_args()
    conservation_pipeline.main(args.config, args.n_cores)



