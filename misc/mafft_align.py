import argparse
import subprocess

def run_mafft(input_file, output_file, num_threads=32):
    # MAFFT command with specified number of threads and first sequence as reference
    mafft_cmd = f"mafft --thread {num_threads} --auto --addfragments {input_file} > {output_file}"

    # Run MAFFT using subprocess
    try:
        subprocess.run(mafft_cmd, shell=True, check=True)
        print("MAFFT alignment completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running MAFFT: {e}")

def main():
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Run MAFFT alignment with specified parameters")

    # Add command-line arguments
    parser.add_argument("input_file", help="Input FASTA file")
    parser.add_argument("output_file", help="Output alignment file")
    parser.add_argument("--num_threads", type=int, default=32, help="Number of threads to use (default: 32)")

    # Parse command-line arguments
    args = parser.parse_args()

    # Run MAFFT alignment function with specified arguments
    run_mafft(args.input_file, args.output_file, args.num_threads)

if __name__ == "__main__":
    main()
