import os
import subprocess
import hashlib
from datetime import datetime

def _get_checksum(file_name):
    """
    Get the MD5 checksum of a file.
    :param file_name: File to check.
    :return: MD5 checksum.
    """
    # Code from Stack Overflow:
    # http://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file

    hash_md5 = hashlib.md5()

    with open(file_name, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)

    return hash_md5.hexdigest()

rule make_md5sum_tab:
    input:
        fofn = rules.make_new_fofn.output.new_fofn,
    output:
        tab = "results/reads_filtered/{sample}/fastq.md5.tsv"
    threads: 1,
    resources:
         mem=120,
         hrs=48,
    run:
        md5_df = pd.read_csv(input.fofn,header=None,names=["CLEANED_PATH"])
        md5_df["CLEANED_MD5"] = md5_df["CLEANED_PATH"].apply(_get_checksum)
        md5_df["CLEANED_SIZE"] = md5_df["CLEANED_PATH"].apply(os.path.getsize)
        md5_df.to_csv(output.tab, sep="\t", index=False)

rule overwrite_sample_fastq_files:
    input:
        cleaned_md5_tab = rules.make_md5sum_tab.output.tab,
    output:
        tab="results/overwrite_records/{sample}.overwrite_records.tab.gz"
    threads: 12,
    resources:
        mem=12,
        hrs=48,
    run:
        fofn_df = get_query_fastq(sample_name=wildcards.sample)
        original_fofn_df = fofn_df.reset_index()
        cleaned_df = pd.read_csv(input.cleaned_md5_tab,sep="\t",header=0)
        cleaned_df["CELL"] = cleaned_df["CLEANED_PATH"].apply(lambda x: x.split('/')[-1].split('.fastq_target')[0])
        original_fofn_df = original_fofn_df.rename(columns={"filepath":"ORIGINAL_PATH"})
        original_fofn_df["CELL"] = original_fofn_df["cell_name"].str.replace(".fastq", "", regex=False)
        original_fofn_df.drop(columns=["cell_name"],inplace=True)

        overwrite_df = pd.merge(original_fofn_df, cleaned_df, on="CELL", how="outer")

        print("Getting MD5s")
        overwrite_df["ORIGINAL_MD5"] = overwrite_df["ORIGINAL_PATH"].apply(_get_checksum)
        overwrite_df["ORIGINAL_SIZE"] = overwrite_df["ORIGINAL_PATH"].apply(os.path.getsize)
        overwrite_df["DATE"] = None
        overwrite_df["STATUS"] = None
        overwrite_df["SAMPLE"] = wildcards.sample
        overwrite_df = overwrite_df[["SAMPLE","CELL", "ORIGINAL_PATH", "ORIGINAL_SIZE", "ORIGINAL_MD5", "CLEANED_PATH", "CLEANED_SIZE", "CLEANED_MD5", "DATE", "STATUS"]]
        print("Copying")
        for index, row in overwrite_df.iterrows():
            original_path = row.ORIGINAL_PATH
            original_md5 = row.ORIGINAL_MD5
            original_tmp_path = f"{row.ORIGINAL_PATH}.tmp"
            cleaned_path = row.CLEANED_PATH
            cleaned_md5 = row.CLEANED_MD5
            error_step = None
            # Checking md5sum
            if original_md5 == cleaned_md5:
                overwrite_df["STATUS"] = "Skipped:Identical"
                overwrite_df["DATE"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                continue
            # Checking file existance
            if (overwrite_df['ORIGINAL_PATH'].isnull() | overwrite_df['CLEANED_PATH'].isnull()).any():
                overwrite_df["STATUS"] = "Skipped:Missing"
                overwrite_df["DATE"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                continue
            # Checking folder permission
            original_path_dir = "/".join(original_path.split('/')[:-1])
            if not os.access(original_path, os.W_OK): # write permission is not granted.
                overwrite_df["STATUS"] = "PermissionDenied"
                overwrite_df["DATE"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                continue

            ## Data overwritting rename (tmp) -> hardlink (or copy) -> re-indexing -> remove tmp 
            
            # Checking temporary file from last overwritting.
            if os.path.isfile(original_tmp_path): # If the original file is incomplete due to an interruption during copying.
                subprocess.run(["mv", original_tmp_path, original_path]) # Recover it from the temporary file which is same as the original.
            # rename
            subprocess.run(["mv", original_path, original_tmp_path])

            # hardlink or copy
            error_flag = 0
            try:
                subprocess.run(["ln", cleaned_path, original_path], check=True, stderr=subprocess.PIPE)
                print ("Generated hard-link")
            except subprocess.CalledProcessError:
                try:
                    subprocess.run(["cp", cleaned_path, original_path], check=True)
                except subprocess.CalledProcessError:
                    error_flag = 1
                    error_step = "datacopy"
                    print (f"###Error found in {error_step} step.")

            # check file migration
            if error_flag == 0: # copy passed
                if os.path.getsize(original_path) == os.path.getsize(cleaned_path): # breifly check the file size after file gets moved.
                    # Re-indexing
                    try:
                        subprocess.run(["bgzip", "-t", original_path], check=True)
                        subprocess.run(["samtools", "faidx", original_path], check=True)
                        subprocess.run(["bgzip", "-r", original_path], check=True)
                    except subprocess.CalledProcessError:
                        error_flag = 1
                        error_step = "reindexing"
                        print (f"###Error found in {error_step} step.")
                    # remove tmp
                    subprocess.run(["rm", "-f", original_tmp_path])
            else: # not successfully moved.
                # recover original file
                subprocess.run(["mv", original_tmp_path, original_path])
                
            # error flag check
            if error_flag == 0: # succeed all process.
                overwrite_df["STATUS"] = "Overwritten"
            else:
                overwrite_df["STATUS"] = "Failed(%s)"%error_step
            overwrite_df["DATE"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        overwrite_df.to_csv(output.tab, sep="\t", index=False, compression="gzip")
        
