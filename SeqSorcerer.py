import os
import subprocess
import glob
import shutil
import sys

#If samtools fails
#conda install bioconda::samtools

#if cutadapt fails install try this 
#pip install cutadapt

#if trim galore install fails try this
#conda install bioconda::trim-galore 


#if feature counts fails then try this: Note that subread contains the featurecounts package
#conda install -c bioconda subread


def check_and_install_dependency(dep_cmd, conda_pkg=None, alt_channel=None):
    if shutil.which(dep_cmd) is None:
        pkg = conda_pkg if conda_pkg else dep_cmd
        install_cmd = ["conda", "install", "-y"]
        if alt_channel:
            install_cmd.extend(["-c", alt_channel])
        install_cmd.append(pkg)
        print(f"[INFO] Dependency '{dep_cmd}' not found. Attempting conda install: {' '.join(install_cmd)}")
        try:
            subprocess.run(install_cmd, check=True)
        except subprocess.CalledProcessError:
            print(f"[ERROR] Failed to install {dep_cmd} via conda. Please install it manually.")
            sys.exit(1)
        if shutil.which(dep_cmd) is None:
            print(f"[ERROR] {dep_cmd} still not found after installation.")
            sys.exit(1)
        print(f"[INFO] Successfully installed {dep_cmd}.")
    else:
        print(f"[INFO] Dependency '{dep_cmd}' found.")

# Check for required dependencies.
# Use the 'bioconda' channel for hisat2 on platforms where it's unavailable on defaults.
dependencies = {
    "hisat2": {"pkg": "hisat2", "alt_channel": "bioconda"},
    "trim_galore": {"pkg": "trim_galore", "alt_channel": None},
    "samtools": {"pkg": "samtools", "alt_channel": "bioconda"},
    "featureCounts": {"pkg": "subread", "alt_channel": "bioconda"}  # featureCounts is part of subread
}
for dep, opts in dependencies.items():
    check_and_install_dependency(dep, opts["pkg"], opts["alt_channel"])

def run_command(command, description):
    cmd_str = " ".join(command) if isinstance(command, list) else command
    print(f"\n[INFO] {description}")
    print(f"[COMMAND] Running: {cmd_str}")
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        output = result.stdout.decode().strip()
        if output:
            print(f"[OUTPUT] {output}")
        return output
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr.decode().strip()
        print(f"[ERROR] {description}: {error_msg}")
        raise

def check_tool_version(tool_name, command_list):
    print(f"\n[INFO] Checking {tool_name} version...")
    print(f"[COMMAND] Running: {' '.join(command_list)}")
    try:
        result = subprocess.run(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        stdout_str = result.stdout.decode('utf-8').strip()
        print(f"[OUTPUT] {tool_name} version: {stdout_str}")
        return stdout_str
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr.decode('utf-8').strip()
        print(f"[ERROR] {tool_name} version check failed: {error_msg}")
        return None

def trim_and_fastqc(folder_path, output_dir):
    print("\n[INFO] Starting trimming and FastQC...")
    r1_files = sorted(glob.glob(os.path.join(folder_path, '*_R1_001.fastq.gz')))
    r2_files = sorted(glob.glob(os.path.join(folder_path, '*_R2_001.fastq.gz')))

    if len(r1_files) != len(r2_files):
        print("[ERROR] Mismatch between R1 and R2 file counts.")
        return

    for r1, r2 in zip(r1_files, r2_files):
        base_name = os.path.basename(r1).replace('_R1_001.fastq.gz', '')
        sample_out = os.path.join(output_dir, base_name)
        os.makedirs(sample_out, exist_ok=True)
        trim_cmd = [
            "trim_galore", "--fastqc", "--paired", "--length", "20",
            "--clip_R2", "15", "--three_prime_clip_R1", "15",
            "-o", sample_out, r1, r2
        ]
        print(f"[INFO] Trimming sample: {base_name}")
        run_command(trim_cmd, f"trim_galore for {base_name}")

def build_reference_genome(fasta_location, genome_dir):
    os.makedirs(genome_dir, exist_ok=True)
    index_prefix = os.path.join(genome_dir, "genome_index")
    build_cmd = ["hisat2-build", "--large-index", "-p", "16", fasta_location, index_prefix]
    run_command(build_cmd, "hisat2-build")
    return index_prefix

def fastagffcheck(fasta_location, gtf_location):
    print("\n[INFO] Performing FASTA and GTF header check...")
    with open(fasta_location, 'r') as fasta_file:
        fasta_header = fasta_file.readline().strip()
        if not fasta_header.startswith('>'):
            raise ValueError("Invalid FASTA file: header must start with '>'")
        fasta_id = fasta_header[1:].split()[0]
    
    gtf_id = None
    with open(gtf_location, 'r') as gtf_file:
        for line in gtf_file:
            if line.startswith('##sequence-region'):
                parts = line.split()
                if len(parts) >= 2:
                    gtf_id = parts[1]
                    break
    if not gtf_id:
        print("[WARN] No '##sequence-region' found in GTF file; header check skipped.")
    elif fasta_id not in gtf_id:
        print(f"[WARN] FASTA ID '{fasta_id}' not found in GTF header '{gtf_id}'.")
    else:
        print("[INFO] FASTA and GTF headers match.")

def automated_alignment(input_folder, output_folder, index_prefix):
    print("\n[INFO] Starting alignment with hisat2...")
    # Recursively find trimmed R1 files in all subdirectories.
    r1_files = glob.glob(os.path.join(input_folder, '**', '*_R1_001_val_1.fq.gz'), recursive=True)
    if not r1_files:
        print("[ERROR] No trimmed R1 files found in the specified folder.")
        return
    for r1_path in r1_files:
        # Derive the sample name from the filename.
        sample_name = os.path.basename(r1_path).replace('_R1_001_val_1.fq.gz', '')
        # Look for the corresponding R2 file in the same directory.
        r2_path = os.path.join(os.path.dirname(r1_path), sample_name + '_R2_001_val_2.fq.gz')
        if not os.path.exists(r2_path):
            print(f"[ERROR] R2 file for sample '{sample_name}' not found in {os.path.dirname(r1_path)}.")
            continue
        summary_file = os.path.join(output_folder, f"{sample_name}_summary.txt")
        bam_file = os.path.join(output_folder, f"{sample_name}.bam")
        hisat2_cmd = f"hisat2 -t --rna-strandness RF --summary-file {summary_file} -p 24 -x {index_prefix} -1 {r1_path} -2 {r2_path}"
        samtools_view_cmd = ["samtools", "view", "-@","24", "-Shu", "-"]
        samtools_sort_cmd = ["samtools", "sort", "-@","24", "-n", "-o", bam_file]
        print(f"\n[INFO] Aligning sample: {sample_name}")
        print(f"[COMMAND] Running alignment pipeline: {hisat2_cmd} | {' '.join(samtools_view_cmd)} | {' '.join(samtools_sort_cmd)}")
        try:
            hisat2_proc = subprocess.Popen(hisat2_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            view_proc = subprocess.Popen(samtools_view_cmd, stdin=hisat2_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            hisat2_proc.stdout.close()
            sort_proc = subprocess.Popen(samtools_sort_cmd, stdin=view_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            view_proc.stdout.close()
            sort_output, sort_error = sort_proc.communicate()
            if sort_proc.returncode != 0:
                print(f"[ERROR] Sorting failed for sample {sample_name}: {sort_error.decode().strip()}")
            else:
                print(f"[INFO] Sample '{sample_name}' aligned successfully.")
        except Exception as e:
            print(f"[ERROR] Alignment failed for sample '{sample_name}': {e}")


def run_feature_counts(gtf_file, input_folder, output_folder):
    print("\n[INFO] Starting gene quantification with featureCounts...")
    for filename in os.listdir(input_folder):
        if filename.endswith('.bam'):
            sample_name = filename.replace('.bam', '')
            bam_file = os.path.join(input_folder, filename)
            output_file = os.path.join(output_folder, f"{sample_name}_counts.txt")
            command = [
                'featureCounts', '-a', gtf_file, '-o', output_file,
                '-T', '2', '-s', '2', '-Q', '0', '-t', 'CDS', '-g', 'transcript_id',
                '--minOverlap', '1', '--fracOverlap', '0', '--fracOverlapFeature', '0',
                '-p', '-C', bam_file
            ]
            print(f"\n[INFO] Running featureCounts for sample: {sample_name}")
            print(f"[COMMAND] Running: {' '.join(command)}")
            try:
                subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
                print(f"[INFO] featureCounts completed for sample '{sample_name}'.")
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] featureCounts failed for sample '{sample_name}': {e.stderr.decode().strip()}")

def rna_seq_pipeline():
    print("\n=== Welcome to SeqSorcerer: Automated RNA-seq Pipeline ===")
    raw_folder = input("Enter the folder path containing your raw FASTQ files (formatted suffix _R1_001.fastq.gz and _R2_001.fastq.gz): ").strip()
    output_dir = input("Enter the desired output directory: ").strip()
    
    trimmed_dir = os.path.join(output_dir, "trimmed_folder")
    os.makedirs(trimmed_dir, exist_ok=True)
    
    aligned_dir = os.path.join(output_dir, "aligned_folder")
    os.makedirs(aligned_dir, exist_ok=True)
    
    counts_dir = os.path.join(output_dir, "counts_folder")
    os.makedirs(counts_dir, exist_ok=True)
    
    gtf_file = input("Enter the location of your reference GTF file: ").strip()
    fasta_file = input("Enter the location of your reference FASTA file: ").strip()

    # Check if input is a GFF file and convert to GTF if needed
    if gtf_file.endswith('.gff'):
        gtf_converted = gtf_file.replace('.gff', '.gtf')
        run_command(["gffread", gtf_file, "-T", "-o", gtf_converted], "Converting GFF to GTF")
        gtf_file = gtf_converted  # Update to use the converted file

    fastagffcheck(fasta_file, gtf_file)

    print("\n[INFO] Indexing reference genome...")
    genome_dir = os.path.join(output_dir, "genome_index")
    index_prefix = build_reference_genome(fasta_file, genome_dir)
    
    print("\n[INFO] Starting trimming and quality control...")
    trim_and_fastqc(raw_folder, trimmed_dir)
    
    print("\n[INFO] Starting alignment...")
    automated_alignment(trimmed_dir, aligned_dir, index_prefix)
    
    print("\n[INFO] Starting gene quantification...")
    run_feature_counts(gtf_file, aligned_dir, counts_dir)
    
    print("\n=== RNA-seq pipeline complete. ===")


if __name__ == "__main__":
    tools = {
        "hisat2": ["hisat2", "--version"],
        "trim_galore": ["trim_galore", "--version"],
        "featureCounts": ["featureCounts", "-v"]
    }
    for tool, cmd in tools.items():
        check_tool_version(tool, cmd)
    
    rna_seq_pipeline()

