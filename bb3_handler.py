from pathlib import Path
from Bio import SeqIO
import xml.etree.ElementTree as ET
import csv
import json

def convert_to_clustal(input_file, output_file):
    """
    Converts a multiple sequence alignment file to clustal format.

    input_file: str, path to the input file
    output_file: str, path to the output file
    """
    SeqIO.convert(input_file, "msf", output_file, "clustal")


#converts all msf files in the input_files directory to clustal format and saves them in the output_files directory
def file_convert_msf_to_clustal(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    for file_path in input_dir.glob("*.msf"):
        output_file = output_dir / (file_path.stem + ".clustal")

        convert_to_clustal(file_path, output_file)

        print(f"Converted {file_path} to {output_file}")

        
    print("All files have been converted to clustal format.")


#file_convert_msf_to_clustal(r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\RV50', r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\RV50\clustal_files')

def parse_bio_xml(file_path):
    """Extracts organism data from a single XML file."""
    tree = ET.parse(file_path)
    root = tree.getroot()
    
    results = []
    for seq in root.findall('.//sequence'):
        info = seq.find('seq-info')
        if info is not None and info.find('organism') is not None:
            entry = {
                "name": seq.find('seq-name').text,
                "organism": info.find('organism').text,
                "definition": info.find('definition').text if info.find('definition') is not None else ""
            }
            results.append(entry)
    return results

def process_folder(folder_path, output_csv):
    """Processes all XMLs and saves to a hybrid CSV."""
    path = Path(folder_path)
    
    with open(output_csv, 'a', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        # Header: File name acts as the index, Data contains the JSON blob
        writer.writerow(["file_id", "protein_data_json"])
        
        for xml_file in path.glob('*.xml'):
            print(f"Processing {xml_file.name}...")
            
            # 1. Parse the data
            parsed_data = parse_bio_xml(xml_file)
            
            # 2. Convert the list/dictionary to a JSON string
            json_blob = json.dumps(parsed_data)
            
            # 3. Write row: index by filename (minus extension)
            writer.writerow([xml_file.stem, json_blob])

process_folder(r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\RV11', r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\bioinformatics_master_list.csv')
process_folder(r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\RV12', r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\bioinformatics_master_list.csv')
process_folder(r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\RV20', r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\bioinformatics_master_list.csv')
process_folder(r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\RV30', r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\bioinformatics_master_list.csv')
process_folder(r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\RV40', r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\bioinformatics_master_list.csv')
process_folder(r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\RV50', r'C:\Users\Sam\Phylogenetic-Tree-Building-Analysis\bb3_release\bioinformatics_master_list.csv')