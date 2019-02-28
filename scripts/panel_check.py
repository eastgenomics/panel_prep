import os
import pwd
import glob
import gzip
import shutil
import subprocess
import MySQLdb

import time
import datetime

import begin
import pysam

import xlrd
import xlwt
from xlutils.copy import copy

style = xlwt.XFStyle()
style.alignment.wrap = 1

genes2transcripts_path = "/mnt/storage/data/NGS/nirvana_genes2transcripts"
genes2transcripts_log = "/mnt/storage/data/NGS/nirvana_genes2transcripts.log"
timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d%H%M%S')


def get_nirvana_dir():
    command = "which nirvana.py"
    wrapper_path = subprocess.check_output(command, shell=True)
    nirvana_path = wrapper_path.split("bin")[0]
    return nirvana_path

def get_nirvana_gff():
    nirvana_dir = get_nirvana_dir()
    nirvana_gff = os.path.join(nirvana_dir, "Data/Cache/24/GRCh37/GRCh37_RefSeq_24.gff.gz")  # May need tweaking for future versions (which don't use 24)
    return nirvana_gff

def check_nirvana_gene_transcript(gene_name, transcript):
    
    gene_name = gene_name.upper()
    nirvana_gff = get_nirvana_gff()
        
    transcripts = set()

    with gzip.open(nirvana_gff) as nir_fh:
        for line in nir_fh:
            gff_gene_name = None
            gff_transcript = None
            gff_tag = ""

        
            fields = line.strip().split("\t")
            record_type = fields[2]
            if not record_type == "transcript":
                continue
            
            info_field = fields[8]
            info_fields = info_field.split("; ")
            for field in info_fields:
                key, value = field.split(" ")
            
                if key == "transcript_id" and value.startswith('"NM_'):
                    gff_transcript = value.replace('"','')

                if key == "gene_name":
                    gff_gene_name = value.replace('"','').upper()
                   
            # This record matches the requested transcript
            if (transcript.strip() == gff_transcript) and (gene_name.strip() == gff_gene_name):
                chrom = fields[0].replace("chr", "")
                start, end = fields[3:5]
                return True, ""

    available_transcripts = "\n".join(nirvana_transcripts(gene_name))
    return False, "Transcript unavailable.\nAvailable transcripts:%s" % available_transcripts

def open_ga_db():
    db = MySQLdb.connect(host="sql01",        # your host, usually localhost
                     user="ccbg_admin",           # your username
                     passwd="ccbg",             # your password
                     db="genetics_ark_1_1_0")   # name of the data base

    return db

def database_query(db, sql):
    # you must create a Cursor object. It will let
    #  you execute all the queries you need
    cur = db.cursor()

    # Use all the SQL you like
    cur.execute(sql)

    # print all the first cell of all the rows
    result = cur.fetchall()
    return result

def close_db(db):
    db.close()

def get_regions(transcript, db):
    regions = {}
    sql = """   SELECT tr.region_id, tr.exon_nr, r.chrom, r.start, r.end, (r.end - r.start + 1) as len 
                FROM transcript AS t 
                LEFT JOIN transcript_region AS tr ON t.id = tr.transcript_id 
                LEFT JOIN region AS r ON tr.region_id = r.id 
                WHERE t.refseq = "{transcript}" ORDER BY tr.exon_nr""".format(transcript=transcript)
    
    #sql = """SELECT r.rid, r2g.exon_nr, r.chr, r.start, r.end, (r.end - r.start +1) as len 
    #            FROM region as r 
    #            INNER JOIN region2gene as r2g on r.rid = r2g.rid 
    #            INNER JOIN gene as g on g.gid = r2g.gid 
    #            WHERE g.name = \"%s\" order by r2g.exon_nr""" % transcript

    region_data = database_query(db, sql)
    
    for record in region_data:
        rid, exon, chrom, start, end, length = record
        if str(chrom).isdigit():
            chrom = int(chrom)

        regions[str(rid)] = {   "exon":int(exon),
                                "chr":chrom,
                                "start":int(start),
                                "end":int(end),
                                "length":int(length)}
    return regions

def combine_depths(run_folder, output_filepath):

    pattern = os.path.join(run_folder, "stats/X*.nirvana_203_5bp.gz")
    depth_files = glob.glob(pattern)

    depth_dict = {}
    
    for depth_file in depth_files:
        
        with gzip.open(depth_file) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue

                chrom, start, end, min_, mean, max_ = line.strip().split()[:6]
                exon_id = "{chrom}-{start}-{end}".format(chrom=chrom, start=start, end=end)
                depth_dict.setdefault(exon_id, {        "min":  0,
                                                        "mean": 0,
                                                        "max":  0,
                                                        "count":0  })
                # Add data from each file to dict
                #pprint.pprint(depth_dict)
                depth_dict[exon_id]["count"] += 1
                depth_dict[exon_id]["min"] += int(min_)
                depth_dict[exon_id]["mean"] += float(mean)
                depth_dict[exon_id]["max"] += int(max_)

    uncompressed_output_filepath = output_filepath.split(".gz")[0]
    with open(uncompressed_output_filepath, "w") as out_fh:
        for exon_id in sorted(depth_dict):
            entry_count = float(depth_dict[exon_id]["count"])
            
            # Calculate mean values for each field across all samples in run
            for field in ["min","mean","max"]:
                depth_dict[exon_id][field]  = round(depth_dict[exon_id][field] / entry_count , 2)
            
            #gene, enst, nm_prefix, refseq, exon_num = exon_id.split("_")
            chrom, start, end = exon_id.split("-")
            min_ = depth_dict[exon_id]["min"]
            mean = depth_dict[exon_id]["mean"]
            max_ = depth_dict[exon_id]["max"]


            output_line = "\t".join(map(str, [chrom, start, end, min_, mean, max_, depth_dict[exon_id]["count"] ]))+"\n"
            out_fh.write(output_line)

    # sort, bgzip and index
    command = "sort -k1V -k2n -o {uncompressed_output_filepath} {uncompressed_output_filepath}".format(uncompressed_output_filepath=uncompressed_output_filepath, output_filepath=output_filepath)
    subprocess.check_call(command, shell=True)
    
    command = "bgzip {uncompressed_output_filepath}".format(uncompressed_output_filepath=uncompressed_output_filepath)
    subprocess.check_call(command, shell=True)
    
    assert os.path.exists(pysam.tabix_index(output_filepath, seq_col=0, start_col=1, end_col=2)), "Unable to generate tabix index"
    
    return output_filepath

def get_coverage_stats(transcript, run_folder, fancy=True, verbose=False):
    coverage = {}

    db = open_ga_db()
    regions = get_regions(transcript, db)
    
    if run_folder.endswith("/"):
        rf_name = os.path.basename(os.path.dirname(run_folder))
    else:
        rf_name = os.path.basename(run_folder)

    rf_depth_file = os.path.join(run_folder, "stats/%s.nirvana_203_5bp.gz" % rf_name)
    
    if not os.path.exists(rf_depth_file):
        combine_depths(run_folder, rf_depth_file)

    subopt_exons = []
    pass_exons = []
    exon_count = 0
    message = ""

    if verbose:
        print " ".join(["chr", "start", "end", "min", "mean", "max", "count"])

    with pysam.TabixFile(rf_depth_file) as tbx:
        for region_id, region_info in regions.items():
            for row in tbx.fetch(region_info["chr"],region_info["start"],region_info["end"]):
                row_chrom, row_start, row_end, row_min, row_mean, row_max, row_count = row.strip().split("\t")
                spaced_row = " ".join(row.strip().split("\t"))

                if  str(row_chrom) == str(region_info["chr"])\
                and int(row_start) == int(region_info["start"])\
                and int(row_end) == int(region_info["end"]):
                    exon_count += 1
                    if float(row_min) < 20:
                        if verbose:
                            if fancy:
                                print hilite(spaced_row, "red")
                            else:
                                print row
                        subopt_exons.append(spaced_row)
                    else:
                        if verbose:
                            print " ".join(spaced_row.split("\t"))
                        pass_exons.append(spaced_row)
    
    assert exon_count, "No exons found for {transcript}".format(transcript=transcript)

    if not subopt_exons:
        if fancy:
            print hilite("Pass", "green", bold=True)
        else:
            print "Pass"
        return True, message

    else:
        if fancy:
            print hilite("Fail", "red", bold=True)
        else:
            print "Fail"
        
        message += "Exons < 20x ({subopt}/{total})\n".format(subopt=len(subopt_exons), total=exon_count)
        message += " ".join(["chr", "start", "end", "min", "mean", "max", "count"]) + "\n"
        
        for exon in subopt_exons:
            if fancy:
                message += hilite(exon, "red") + "\n"
            else:
                message += exon + "\n"
        
        percent_20x = get_20x_percentage(transcript, run_folder)
        message = "Coverage >=20X: %.2f%%\r\n" % (percent_20x*100) + message
        print message
        return False, message

def transcript2regions(transcript):
    sql_command =  "select id, chrom, start, end from region where id in (select region_id from transcript_region where transcript_id in (select id from transcript where refseq = '%s'))" % transcript
    db = MySQLdb.connect(   host="sql01",    # your host, usually localhost
                            user="ccbg_admin",         # your username
                            passwd="ccbg",  # your password
                            db="genetics_ark_1_1_0")        # name of the data base

    cur = db.cursor()

    # Use all the SQL you like
    cur.execute(sql_command)

    # print all the first cell of all the rows
    regions = []
    for row in cur.fetchall():
        region = row[0:4]
        # Make coords ints
        region = [int(row[0])] + [row[1]] + map(int, row[2:4])
        regions.append(region)
  
    db.close()
    return regions

def get_region_coverage(chrom, start, end, depths_file):
    coverage = {"chrom":                None,
                "start":                None,
                "end":                  None,
                "length":               None,
                "min":                  None,
                "mean":                 None,
                "max":                  None,
                "depth_0_regions":      None,
                "depth_1_5_regions":    None,
                "depth_6_9_regions":    None,
                "depth_10_19_regions":  None,
                "depth_20_regions":     None,
                "depth_0_length":       None,
                "depth_1_5_length":     None,
                "depth_6_9_length":     None,
                "depth_10_19_length":   None,
                "depth_20_length":      None,
                "stdev":                None
                }
    
    # Get all records intersecting our region
    command = "tabix %s %s:%s-%s" % (depths_file, chrom, start, end)
    #logger.debug(command)
    coverage_records = subprocess.check_output(command, shell=True).strip().split("\n")
    
    # Find the region with the same start + end
    for coverage_record in coverage_records:
        record_fields = coverage_record.strip().split("\t")
        #logger.debug("Record fields: %s", ", ".join(record_fields))
        record_start, record_end = map(int, record_fields[1:3])
        
        # Mismatched coords are skipped
        if (start != record_start) or (end != record_end):
            #logger.debug("Region range mismatch\nStart: %s/%s\nEnd: %s/%s\nSkipping this record" % (start, record_start, end, record_end))
            continue
        
        # Different fields for sample vs run folder depths
        if ".nirvana_203_5bp." in os.path.basename(depths_file):  # Sample
            depth_file_type = "sample"
            field_names = "chrom", "start", "end", "min", "mean", "max", "depth_0_regions", "depth_1_5_regions", "depth_6_9_regions", "depth_10_19_regions", "depth_20_regions"    
        
        elif ".refseq_nirvana_5bp." in os.path.basename(depths_file):  # Run
            depth_file_type = "run"
            field_names = ["chrom", "start", "end", "mean", "stdev"]

        while len(record_fields) < len(field_names):  # Trailing empty fields trimmed by tabix - restore these
            record_fields.append("")

        for index, field_value in enumerate(record_fields):
            #logger.debug([index, field_names[index], field_value])
            coverage[field_names[index]] = field_value

        coverage["start"] = int(coverage["start"]) - 5
        coverage["end"] = int(coverage["end"]) + 5
        coverage["length"] = (coverage["end"] - coverage["start"]) + 1

        # Calculate total length of suboptimal depths by depth bin
        if depth_file_type == "sample":
            for depth_bin in "depth_0", "depth_1_5", "depth_6_9", "depth_10_19", "depth_20":
                low_depth_region_length_for_bin = 0
                
                low_depth_regions = coverage[depth_bin+"_regions"]
                if low_depth_regions:
                    low_depth_regions = coverage[depth_bin+"_regions"].split(",")
                
                # Sum the lengths of all sub-regions in that bin
                for low_depth_region in low_depth_regions:
                    low_depth_region_start, low_depth_region_end = map(int, low_depth_region.split(":")[-1].split("-"))  # 1:1234-5678 ---> 1234, 5678
                    low_depth_region_length = (low_depth_region_end - low_depth_region_start) + 1
                    low_depth_region_length_for_bin += low_depth_region_length

                coverage[depth_bin+"_length"] = low_depth_region_length_for_bin

        return coverage

def get_20x_percentage(transcript, run_folder, gene=None, verbose=False):
    if verbose:
        sys.stderr.write('%s\n' % gene)

    run_folder = os.path.basename(run_folder)
    run_depth_file = "/mnt/storage/data/NGS/%s/stats/%s.refseq_nirvana_5bp.gz" % (run_folder, run_folder)
    sample_depth_files = glob.glob("/mnt/storage/data/NGS/%s/stats/*nirvana_203_5bp.gz" % run_folder)[0:20]
    genes2transcripts_filepath = "/mnt/storage/data/NGS/nirvana_genes2transcripts"

    regions = transcript2regions(transcript)
    
    gene_coverage = {"gene":                 gene,  #Done
                     "transcript":           transcript,
                     "chrom":                None,  #Done
                     "start":                float('inf'),    #Done
                     "end":                  -1,      #Done
                     "total_length":         0,     #Done
                     "coding_length":        0,     #Done
                     "min":                  None,  
                     "mean":                 None,
                     "max":                  None,
                     #"depth_0_regions":      None,
                     #"depth_1_5_regions":    None,
                     #"depth_6_9_regions":    None,
                     #"depth_10_19_regions":  None,
                     #"depth_20_regions":     None,
                     "depth_0_length":       None,
                     "depth_1_5_length":     None,
                     "depth_6_9_length":     None,
                     "depth_10_19_length":   None,
                     "depth_20_length":      0,
                     "perc_20x":             None,
                     #"stdev":                None,
                     "regions":              {},
                    }
    
    for region in regions:
        rid, chrom, start, end = region
        #gene_coverage.set_default("chrom", chrom)
        
        #if start <= gene_coverage["start"]:
        #    gene_coverage["start"] = start

        coverage = get_region_coverage(chrom, start, end, run_depth_file)
        gene_coverage["regions"][rid] = coverage

        sample_coverage_dicts = [get_region_coverage(chrom, start, end, sample_depth_file) for sample_depth_file in sample_depth_files]
        number_of_samples = len(sample_coverage_dicts)

        start = int(start) - 5
        end = int(end) + 5
        total_bases_over_20x = (end - start) + 1
        for depth_bin in "depth_0", "depth_1_5", "depth_6_9", "depth_10_19":
            total_suboptimal_bases = sum([sample_coverage_dict[depth_bin+"_length"] for sample_coverage_dict in sample_coverage_dicts])
            suboptimal_bases_per_sample = total_suboptimal_bases/float(number_of_samples)

            #print depth_bin, total_suboptimal_bases, suboptimal_bases_per_sample
            total_bases_over_20x -= suboptimal_bases_per_sample
        
        gene_coverage["depth_20_length"] += total_bases_over_20x
        #print total_bases_over_20x
        #print gene_coverage["regions"][rid]["length"]
        #print chrom, start, end

    gene_start = int(min([gene_coverage["regions"][rid]["start"] for rid in gene_coverage["regions"].keys()]))
    #print gene_start
    gene_coverage["start"] = gene_start

    gene_end = int(max([gene_coverage["regions"][rid]["end"] for rid in gene_coverage["regions"].keys()]))
    #print gene_end
    gene_coverage["end"] = gene_end

    gene_chrom = list(set([gene_coverage["regions"][rid]["chrom"] for rid in gene_coverage["regions"].keys()]))
    
    gene_chrom = gene_chrom[0]
    #print gene_chrom
    gene_coverage["chrom"] = gene_chrom

    total_gene_length = (gene_end - gene_start) + 1
    #print total_gene_length
    gene_coverage["total_length"] = total_gene_length

    coding_gene_length = sum([(int(gene_coverage["regions"][rid]["end"]) - int(gene_coverage["regions"][rid]["start"])) + 1 for rid in gene_coverage["regions"].keys()])
    #print coding_gene_length
    gene_coverage["coding_length"] = coding_gene_length    

    min_depth = min([gene_coverage["regions"][rid]["min"] for rid in gene_coverage["regions"].keys()])
    #print gene_coverage["regions"][rid].items()    
        
    mean_depth = sum([gene_coverage["regions"][rid]["length"] * float(gene_coverage["regions"][rid]["mean"]) for rid in gene_coverage["regions"].keys()])/gene_coverage["coding_length"]
    #print mean_depth
    gene_coverage["mean"] = mean_depth

    percentage_20x = gene_coverage["depth_20_length"]/gene_coverage["coding_length"]
    gene_coverage["perc_20x"] = percentage_20x

    if verbose:
        print "\t".join(map(str, [gene_coverage["gene"], gene_coverage["transcript"], gene_coverage["coding_length"], gene_coverage["mean"], gene_coverage["perc_20x"]]))

    return gene_coverage["perc_20x"]

def get_panels_containing_gene(gene):
    genepanels = "/mnt/storage/data/NGS/genepanels"
    panels_containing_gene = []

    with open(genepanels) as gp_fh:
        for line in gp_fh:
            if line.startswith("#"):
                continue
            
            fields = line.strip().split("\t")
            if len(fields) != 3:
                continue
            panel, pid, gene_symbol = fields
            if gene == gene_symbol:
                panels_containing_gene.append(panel)
    return panels_containing_gene

def hilite(string, colour, bold=False):
    """
    Colour tty text output
    
    Args:
        string (str): The string to be highlighted
        colour (str): A color. "red" or "green"
        bold (bool, optional): True if bold formatting required (default=False)
    
    Returns:
        str: String which appears highlighed as specified when printed in tty
    """
    
    attr = []
    if colour == "green":
        attr.append('32')
    elif colour == "red":
        attr.append('31')
    if bold:
        attr.append('1')
    return '\x1b[%sm%s\x1b[0m' % (';'.join(attr), string)

def get_username():
    return pwd.getpwuid( os.getuid() )[ 0 ]

@begin.subcommand
def archive_genes2transcripts(genes2transcripts_path=genes2transcripts_path):
    archive_dir = "/mnt/storage/data/NGS/nirvana_genes2transcripts_archive"
    archive_filename = os.path.basename(genes2transcripts_path) + "_" + timestamp
    archive_file_path = os.path.join(archive_dir, archive_filename)
    shutil.copyfile(genes2transcripts_path, archive_file_path)
    return archive_file_path

@begin.subcommand
def panel_definition_file(form):

    form_dir = os.path.dirname(form)
    output_filepath = os.path.join(form_dir, "NewLargePanel.txt")

    panel_group_ids = { "Blood":1,
                        "Cancer":2,
                        "Cardiac":3,
                        "Dermatology":4,
                        "Developmental Delay":5,
                        "Dysmorphology":6,
                        "Eye Disorders":7,
                        "Gastrointestinal":8,
                        "Hearing":9,
                        "Immunology":10,
                        "Liver":11,
                        "Metabolic":12,
                        "Neurology":13,
                        "Neuromuscular":14,
                        "Pancreatic":15,
                        "Pulmonary":16,
                        "Renal":17,
                        "Sexual Development":18,
                        "Skeletal":19,
                        "Whole Exome":20,
                        "BLANK":21,
                        "Neurometabolic":22,
                        "Chromosome Instability":23,
                        "Nuclear Mitochondrial":24,
                        }
    rb = xlrd.open_workbook(form)
    rs = rb.sheet_by_name("Sheet1")
    panel_name = rs.cell(2,2).value
    panel_group_name = rs.cell(1,2).value

    panel_group_id = panel_group_ids.get(panel_group_name, None)

    assert panel_group_id, "No panel group id for '{panel_name}'".format(panel_name=panel_name)
    
    header_row_index = 8

    header_row = rs.row(header_row_index)

    header_row_values = [cell.value for cell in header_row]

    gene_symbol_column_index    = header_row_values.index("Gene symbol")
    hgnc_id_column_index        = header_row_values.index("HGNC ID")
    transcript_column_index     = header_row_values.index("Transcript")
    symboltx_check_column_index   = header_row_values.index("Symbol / Transcript check")
    coverage_check_column_index = header_row_values.index("20X Coverage check")
    g2t_check_column_index      = header_row_values.index("G2T check")
    bix_notes_column_index      = header_row_values.index("Bioinformatics Notes")

    nrows = rs.nrows

    genes=set()

    for row_index in range(header_row_index+1, nrows-1):
        row = rs.row(row_index)
        row_values = [cell.value for cell in row]
        gene_symbol = row_values[gene_symbol_column_index]
        symboltx_check = row_values[symboltx_check_column_index]
        coverage_check = row_values[coverage_check_column_index]
        g2t_check = row_values[g2t_check_column_index]
        
        if not gene_symbol:
            continue
        
        if all([gene_symbol, symboltx_check, coverage_check, g2t_check]):
            genes.add(gene_symbol)
        else:
            exit("Bailing due to failed checks in %s" % gene)

    with open(output_filepath, "w") as out_fh:
        out_fh.write(panel_name+"\r\n")
        out_fh.write(str(panel_group_id)+"\r\n")
        for gene in sorted(list(genes)):
            out_fh.write(gene+"\r\n")

    print "\nPanel definition file: {output_filepath}".format(output_filepath=output_filepath)

@begin.subcommand
def add_to_genes2transcripts(gene_name, transcript, genes2transcripts_path=genes2transcripts_path, genes2transcripts_log=genes2transcripts_log):
    assert archive_genes2transcripts()

    assert transcript.startswith("NM_"), "Transcript must be Refseq"

    transcript_available, transcript_available_message = check_nirvana_gene_transcript(gene_name, transcript)
    genes2transcripts_check, genes2transcripts_message = check_genes2transcripts(gene_name, transcript, silent=True)
    
    if genes2transcripts_check:
        print "{gene}\t{transcript} already present in {genes2transcripts_path}".format(gene=gene_name, transcript=transcript, genes2transcripts_path=genes2transcripts_path)
        return False

    # It can be added, but it's currently missing. Try to add it.
    if transcript_available and (genes2transcripts_message == "Gene and transcript missing from genes2transcripts"):

        command = 'echo -e "{gene}\t{transcript}" >> {genes2transcripts}'.format(gene=gene_name, transcript=transcript, genes2transcripts=genes2transcripts_path)
        
        try:
            subprocess.check_call(command, shell=True)
            print command
        except subprocess.CalledProcessError as error:
            print error
            return False

        # Record in log
        username  = get_username()
        command   = 'echo -e "{timestamp}: {username} added {gene_name}\t{transcript}" >> {genes2transcripts_log}'.format(timestamp=timestamp, username=username, gene_name=gene_name, transcript=transcript, genes2transcripts_log=genes2transcripts_log)
        if subprocess.check_call(command, shell=True):  # Returns 0 if OK
            print "Unable to log nirvana_genes2transcripts changes"
        return True

@begin.subcommand
def remove_from_genes2transcripts(gene_name, transcript=None, genes2transcripts_path=genes2transcripts_path, genes2transcripts_log=genes2transcripts_log):
    assert archive_genes2transcripts()
    if transcript:
        assert transcript.startswith("NM_"), "Transcript must be Refseq"

    lines = []
    with open(genes2transcripts_path) as g2t_fh:
        for line in g2t_fh:
            fields = line.strip().split("\t")
            if fields[0] == gene_name:
                continue
            lines.append(line)

    with open(genes2transcripts_path, "w") as g2t_fh:
        for line in lines:
            g2t_fh.write(line)

    # Record in log
    username  = get_username()
    command   = 'echo -e "{timestamp}: {username} removed {gene_name}" >> {genes2transcripts_log}'.format(timestamp=timestamp, username=username, gene_name=gene_name, genes2transcripts_log=genes2transcripts_log)
    if subprocess.check_call(command, shell=True):  # Returns 0 if OK
        print "Unable to log nirvana_genes2transcripts changes"
    return True

    
   

    if transcript_available and (genes2transcripts_message == "Gene and transcript missing from genes2transcripts"):

        command = 'echo -e "{gene}\t{transcript}" >> {genes2transcripts}'.format(gene=gene_name, transcript=transcript, genes2transcripts=genes2transcripts_path)
        
        try:
            subprocess.check_call(command, shell=True)
            print command
        except subprocess.CalledProcessError as error:
            print error
            return False

        # Record in log
        username  = get_username()
        command   = 'echo -e "{timestamp}: {username} added {gene_name}\t{transcript}" >> {genes2transcripts_log}'.format(timestamp=timestamp, username=username, gene_name=gene_name, transcript=transcript, genes2transcripts_log=genes2transcripts_log)
        if subprocess.check_call(command, shell=True):  # Returns 0 if OK
            print "Unable to log nirvana_genes2transcripts changes"
        return True

@begin.subcommand
def update_genes2transcripts(gene_name, transcript, genes2transcripts_path=genes2transcripts_path, genes2transcripts_log=genes2transcripts_log):
    remove_from_genes2transcripts(gene_name)
    add_to_genes2transcripts(gene_name, transcript)
    return True

@begin.subcommand
def check_genes2transcripts(gene_symbol, transcript, fancy=False, silent=False):

    missing = True
    result = False

    with open(genes2transcripts_path) as g2t_fh:
        for line in g2t_fh:
            
            if line.startswith("#"):
                continue
            
            fields = line.strip().split("\t")
            if len(fields) != 2:
                continue
            
            message = ""
            g2t_gene, g2t_transcript = fields

            same_symbol = (gene_symbol == g2t_gene)          
            same_transcript = (transcript == g2t_transcript)

            if same_symbol and same_transcript:
                missing = False
                result = True
                if not silent:
                    if fancy:
                        print hilite("Pass", "green", bold=True)
                    else:
                       print "Pass"
                return result, message

            elif same_symbol:  # and diff transcript
                missing = False
                if not silent:
                    if fancy:
                        print hilite("Fail", "red", bold=True)
                    else:
                        print "Fail"
                    
                affected_panels = get_panels_containing_gene(gene_symbol)
                message += "Different transcript already assigned to gene symbol\r\n"
                message += " ".join(line.strip().split("\t"))+"\r\n"
                message += "Changing this transcript will affect the following panels:\r\n"
                message += "\n".join(affected_panels)+"\r\n"
                
                if not silent:
                    print message
                
                return result, message

            elif same_transcript:  # and diff symbol
                missing = False
                if not silent:
                    if fancy:
                        print hilite("Fail", "red", bold=True)
                    else:
                        print "Fail"

                message += "Transcript already assigned to different gene symbol \r\n"
                message += " ".join(fields)
                
                if not silent:
                    print message
                
                return result, message
    
    if missing:
        message = ""
        if not silent:
            if fancy:
                print hilite("Fail", "red", bold=True)
            else:
                print "Fail"
            
        message += "Gene and transcript missing from genes2transcripts"
        
        if not silent:
            print message
        
        return result, message

@begin.subcommand
def transcript_coverage(transcript, run_folder="/mnt/storage/data/NGS/180904_K00178", fancy=True, verbose=False):
    return get_coverage_stats(transcript, run_folder, verbose=verbose, fancy=fancy)

@begin.subcommand
def panel_check(gene_symbol, transcript, run_folder="/mnt/storage/data/NGS/180904_K00178", fancy=True):
    
    assert transcript.startswith("NM_"), "Transcript must be RefSeq (NM_)"
    
    gene_symbol = gene_symbol.strip()
    transcript = transcript.strip()
    
    print
    print gene_symbol, transcript
    
    transcript_available = check_nirvana_gene_transcript(gene_symbol, transcript)
    
    print "Symbol/Tx check:\t", 
    if transcript_available:
        if fancy:
            print hilite("Pass", "green", bold=True)
        else:
            print "Pass"
    else:
        if fancy:
            print hilite("Fail", "red", bold=True)
        else:
            print "Fail"

    print "Coverage check:\t\t", 
    coverage_stats = get_coverage_stats(transcript, run_folder, fancy)
    
    print "g2t check:\t\t",
    g2t_check = check_genes2transcripts(gene_symbol, transcript, fancy)

@begin.subcommand
def check_form(form):

    path, ext = os.path.splitext(form)
    output_filepath = "{path}_{timestamp}{ext}".format(path=path, timestamp=timestamp, ext=ext)

    rb = xlrd.open_workbook(form, formatting_info=1)
    r_sheet = rb.sheet_by_name("Sheet1")
    
    wb = copy(rb)
    w_sheet = wb.get_sheet(0)

    header_row_index = 8

    header_row = r_sheet.row(header_row_index)

    header_row_values = [cell.value for cell in header_row]

    gene_symbol_column_index    = header_row_values.index("Gene symbol")
    hgnc_id_column_index        = header_row_values.index("HGNC ID")
    transcript_column_index     = header_row_values.index("Transcript")
    symboltx_check_column_index   = header_row_values.index("Symbol / Transcript check")
    #transcript_check_column_index = header_row_values.index("Transcript check")
    coverage_check_column_index = header_row_values.index("20X Coverage check")
    g2t_check_column_index      = header_row_values.index("G2T check")
    bix_notes_column_index      = header_row_values.index("Bioinformatics Notes")

    nrows = r_sheet.nrows
    
    result_map = {True:  "Pass",
                  False: "FAIL"}

    for row_index in range(header_row_index+1, nrows-1):
        row = r_sheet.row(row_index)
        row_values = [cell.value for cell in row]

        if any(row_values):
            gene_symbol = row_values[gene_symbol_column_index]
            transcript  = row_values[transcript_column_index]
            
            gene_tx_check_bool, gene_tx_check_message     = check_nirvana_gene_transcript(gene_symbol, transcript)
            
            # If the tx is not in nirvana then we won't have coverage for the region, so don't check for it
            if gene_tx_check_bool:
                coverage_check_bool, coverage_check_message   = transcript_coverage(transcript, fancy=False)
            else:
                coverage_check_bool, coverage_check_message = False, ""
            
            g2t_check_bool, g2t_check_message             = check_genes2transcripts(gene_symbol, transcript)

            gene_tx_check_result = result_map[gene_tx_check_bool]
            coverage_check_result = result_map[coverage_check_bool]
            g2t_check_result = result_map[g2t_check_bool]

            # If missing from g2t and everything else is OK then add it and rerun the check
            if all([gene_tx_check_bool, coverage_check_bool]) and not g2t_check_bool:
                added_to_g2t = add_to_genes2transcripts(gene_symbol, transcript)
                g2t_check_bool, g2t_check_message     = check_genes2transcripts(gene_symbol, transcript)
                g2t_check_result = result_map[g2t_check_bool]

            row_notes = "\n".join([gene_tx_check_message,coverage_check_message,g2t_check_message])
            
            w_sheet.write(row_index, symboltx_check_column_index, gene_tx_check_result)
            w_sheet.write(row_index, coverage_check_column_index, coverage_check_result)
            w_sheet.write(row_index, g2t_check_column_index, g2t_check_result)
            w_sheet.write(row_index, bix_notes_column_index, row_notes, style)

    print output_filepath
    wb.save(output_filepath)

@begin.subcommand
def nirvana_transcripts(gene_name, verbose=True):
    gene_name = gene_name.upper()
    nirvana_dir = get_nirvana_dir()
    nirvana_gff = os.path.join(nirvana_dir, "Data/Cache/24/GRCh37/GRCh37_RefSeq_24.gff.gz")  # May need tweaking for future versions (which don't use 24)
    #nirvana_gff = "/mnt/storage/apps/software/nirvana/2.0.3/Data/Cache/24/GRCh37/GRCh37_RefSeq_24.gff.gz"
    #hgmd = "/data/gemini/HGMD/20170718/HGMD_PRO_2017.2_hg19.vcf"
    transcripts = set()

    with gzip.open(nirvana_gff) as nir_fh:
        for line in nir_fh:
            gff_gene_name = None
            gff_transcript = None
            gff_tag = ""

            if gene_name in line.upper():
                fields = line.strip().split("\t")
                record_type = fields[2]
                if not record_type == "transcript":
                    continue
                info_field = fields[8]
                info_fields = info_field.split("; ")
                for field in info_fields:
                    key, value = field.split(" ")
                
                    if key == "gene_name":
                        gff_gene_name = value.replace('"','').upper()
                    elif key == "transcript_id" and value.startswith('"NM_'):
                        gff_transcript = value.replace('"','')
                    elif key == "tag":
                        gff_tag = value.replace('"','')

                if gff_gene_name == gene_name and gff_transcript:
                    chrom = fields[0].replace("chr", "")
                    start, end = fields[3:5]
                    #print "\t".join([gff_transcript, chrom, start, end, gff_tag])

                    transcripts.add("\t".join([gff_gene_name, gff_transcript, chrom, start, end, gff_tag]))
    
    if verbose:
        for transcript in sorted(transcripts):
            print transcript
    
    return [x.split("\t")[1] for x in transcripts]

@begin.subcommand
def get_HGMD_transcript(gene_name):
    db = open_ga_db()
    sql_request = "SELECT refcore FROM hgmd_pro.gene2refseq WHERE hgmdID = (SELECT gene_id FROM hgmd_pro.allgenes WHERE gene = '{}')".format(gene_name)
    HGMD_transcript = database_query(db, sql_request)
    return HGMD_transcript

@begin.start
def main():
    pass