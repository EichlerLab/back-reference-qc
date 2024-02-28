import os
import re
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

def _file_size_check(original_path, cleaned_path):
    # 1: copied / 0: not copied or interrupted
    if os.path.getsize(original_path) == os.path.getsize(cleaned_path):
        return 1
    else:
        print (f"### SizeDiff After copying: \n Original Path: {original_path}\n Cleaned Path: {cleaned_path}")
        return 0

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

rule make_overwrite_target_list:
    input:
        cleaned_md5_tab = rules.make_md5sum_tab.output.tab,
    output:
        tab="results/overwrite_target_lists/{sample}.overwrite_target_lists.tab.gz"
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
        overwrite_df["SAMPLE"] = wildcards.sample
        overwrite_df = overwrite_df[["SAMPLE","CELL", "ORIGINAL_PATH", "ORIGINAL_SIZE", "ORIGINAL_MD5", "CLEANED_PATH", "CLEANED_SIZE", "CLEANED_MD5"]]
        overwrite_df.to_csv(output.tab, sep="\t", index=False, compression="gzip")
        os.system(f"touch results/overwrite_target_lists/.{wildcards.sample}.overwrite_prepare.done")

rule overwrite_sample_fastq_files:
    input:
        cleaned_flag = "results/overwrite_target_lists/.{sample}.overwrite_prepare.done",
    output:
        tab="results/overwrite_records/{sample}.overwrite_records.tab.gz"
    threads: 12,
    resources:
        mem=12,
        hrs=48,
    run:
        overwrite_df = pd.read_csv(f"results/overwrite_target_lists/{wildcards.sample}.overwrite_target_lists.tab.gz", sep="\t", header=0)
        # default set for DATA / STATUS / INFO
        overwrite_df["DATE"] = None
        overwrite_df["STATUS"] = None
        overwrite_df["INFO"] = None
        print("Copying...")
        for index, row in overwrite_df.iterrows():
            sample_name = overwrite_df.at[index, "SAMPLE"]
            original_path = row.ORIGINAL_PATH
            fastq_base_name = os.path.basename(original_path)
            original_fai = f"{original_path}.fai"
            original_gzi = f"{original_path}.gzi"
            original_md5 = row.ORIGINAL_MD5
            original_tmp_path = f"{row.ORIGINAL_PATH}.tmp"
            cleaned_path = row.CLEANED_PATH
            cleaned_fai = f"{cleaned_path}.fai"
            cleaned_gzi = f"{cleaned_path}.gzi"
            cleaned_md5 = row.CLEANED_MD5
            error_step = None
            permission_check = []
            copy_info = []
            error_flag = 0
            copy_flag = 0

            # Checking file existance
            if pd.isnull(row['ORIGINAL_PATH']) or pd.isnull(row['CLEANED_PATH']):
                    overwrite_df.at[index, "STATUS"] = "Skipped:Missing"
                    overwrite_df.at[index, "DATE"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    continue

            # Checking folder permission
            original_path_dir = "/".join(original_path.split('/')[:-1])

## permission for both of directory and file is needed to move file (generating tmp file).

            if not os.access(original_path_dir, os.W_OK): # not granted for writing in dir.
                permission_check.append("dir")
            if not os.access(original_path, os.W_OK): # not granted for the file.
                permission_check.append("file")
            permission_result = "/".join(permission_check)
            if not permission_result == "": # write permission is not granted.
                overwrite_df.at[index, "STATUS"] = f"Failed(PermissionDenied;{permission_result})"
                overwrite_df.at[index, "DATE"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print (f"### Permission Denied for {permission_result}: {sample_name}-{fastq_base_name}")
                continue

## Data overwritting rename (tmp) -> hardlink (or copy) -> re-indexing -> remove tmp 

            # Checking temporary file from last overwritting.
            if os.path.isfile(original_tmp_path): # If the original file is incomplete due to an interruption during copying.
                subprocess.run(["mv", original_tmp_path, original_path]) # Recover it from the temporary file which is same as the original.
            # Checking md5sum
            if original_md5 == cleaned_md5: # assume copying has completed successfully in advance ( before the creation of target list ).
                copy_flag = _file_size_check(original_path, cleaned_path) # breifly check after copying; 1(copied) or 0(not copied)
                if copy_flag == 0: # MD5 of original fastq has been changed.
                    error_step = "filematch"
                    print (f"###The file {original_path} no longer match to the cleaned fastq.")
                    overwrite_df.at[index, "STATUS"] = f"Failed({error_step})"
                    overwrite_df.at[index, "DATE"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    continue
            else:
                # Secondary check for possibility that the file was already replaced in previous after creating the target list.
                current_original_md5 = _get_checksum(original_path)
                if current_original_md5 == cleaned_md5: # already replaced. Skip without copying.
                    copy_flag = 1
                    print (f"### Already replaced after the creation of target list: {sample_name}-{fastq_base_name}")
                # skipping copy process when it's already done.
                if copy_flag == 0: # not copied yet
                    # rename for tmp fastq
                    try:
                        subprocess.run(["mv", original_path, original_tmp_path])

                    except subprocess.CalledProcessError as e:
                        error_flag = 1
                        error_step = "data_copy"
                        print (f"###Error found in {error_step} step: {sample_name}-{fastq_base_name}")
                    # hardlink or copy
                    try:
                        subprocess.run(["ln", cleaned_path, original_path], check=True, stderr=subprocess.PIPE)
                        copy_flag = _file_size_check(original_path, cleaned_path) # breifly check after copying; 1(copied) or 0(not copied)
                    except:
                        try:
                            subprocess.run(["cp", cleaned_path, original_path], check=True, stderr=subprocess.PIPE)
                            copy_flag = _file_size_check(original_path, cleaned_path) # breifly check after copying; 1(copied) or 0(not copied)
                        except subprocess.CalledProcessError as e:
                            print (f"#CopyError: {e.stderr.decode()}")
                            error_flag = 1
                            error_step = "data_copy"
                            print (f"###Error found in {error_step} step: {sample_name}-{fastq_base_name}")
            if copy_flag == 0: # data copy failed or interrupted.
                error_step = "data_copy"
                # Recover it from the temporary file if it's possible.
                if os.path.isfile(original_tmp_path):
                    subprocess.run(["mv", original_tmp_path, original_path])
                overwrite_df.at[index, "STATUS"] = f"Failed({error_step})"
                overwrite_df.at[index, "DATE"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                continue
## Assuming data copy has been completed successfully.
## Index Copy
            index_pairs = {
                cleaned_fai:original_fai,
                cleaned_gzi:original_gzi,
            }
            try:
                for cleaned_index in index_pairs:
                    cleaned_fastq = '.'.join(cleaned_index.split('.')[:-1])
                    original_index = index_pairs[cleaned_index]
                    print (cleaned_fastq)
                    print (original_index)



                    if os.path.isfile(original_index): # if index exists; need to delete for preventing hardlink affection.
                        if not os.access(original_index, os.W_OK): # write permission denied
                            error_flag = 1
                            error_step = "Permission_Denied(index)"
                            print (f"### Permission Denied for indexing: {sample_name}-{fastq_base_name}")
                            continue
                        else: # write permisson granted
                            try:
                                subprocess.run(["rm", "-f", original_index], check=True, stderr=subprocess.PIPE)
                            except:
                                error_flag = 1
                                error_step = "index_cleaning"
                                print (f"###Error found in {error_step} step: {sample_name}-{fastq_base_name}")
                                continue

## copying indexing after clearing index.
                    # case 1: cleaned index exists.
                    if os.path.isfile(cleaned_index):
                        subprocess.run(["cp", cleaned_index, original_index], check=True, stderr=subprocess.PIPE) # folder and file permisson was already granted.
                    # case 2: cleaned index doesn't exists; need to re-generate for cleaned fastq
                    else:
                        # re-index for cleaned fastq.
                        if re.search('\.gzi$', cleaned_index): # .gzi index
                            try:
                                subprocess.run(["samtools", "fqidx", cleaned_fastq], check=True, stderr=subprocess.PIPE)
                            except:
                                error_flag = 1
                                error_step = "gzi_indexing"
                                print (f"###Error found in {error_step} step: {sample_name}-{fastq_base_name}")
                                continue

                        elif re.search('\.fai$', cleaned_index): # .fai index
                            try:
                                subprocess.run(["samtools", "faidx", cleaned_fastq], check=True, stderr=subprocess.PIPE)
                            except:
                                error_flag = 1
                                error_step = "fai_indexing"
                                print (f"###Error found in {error_step} step: {sample_name}-{fastq_base_name}")
                                continue
                        subprocess.run(["cp", cleaned_index, original_index], check=True, stderr=subprocess.PIPE) # folder permisson was already granted.
                        
            except: # Unexpected error
                error_flag = 1
                error_step = "index_copy(unknown)"
                print (f"###Error found in {error_step} step: {sample_name}-{fastq_base_name}")
## Final summary step.
            if error_flag == 0: # succeed all process without any error.
                # remove tmp fastq
                if os.path.isfile(original_tmp_path):
                    subprocess.run(["rm", "-f", original_tmp_path])
                overwrite_df.at[index, "STATUS"] = "Overwritten"
            else:
                overwrite_df.at[index, "STATUS"] = f"Failed({error_step})"
            overwrite_df.at[index, "DATE"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        overwrite_df.to_csv(output.tab, sep="\t", index=False, compression="gzip")
