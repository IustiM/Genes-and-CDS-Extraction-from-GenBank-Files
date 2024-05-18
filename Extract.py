#M.I.
from Bio import SeqIO

def extract_genes_and_cds(input_file, output_file):
    header = []
    features_section = False

    # Deschidem fisierul de intrare pentru citire
    with open(input_file, "r") as gb_file:
        # Citim si stocam header-ul pâna la sectiunea FEATURES
        for line in gb_file:
            header.append(line)
            if line.startswith("FEATURES"):
                features_section = True
                break

        # Reseteaza pozitia pointerului la inceputul fisierului pentru a citi din nou
        gb_file.seek(0)
        # Citim fisierul GenBank
        record = SeqIO.read(gb_file, "genbank")

    # Deschidem fisierul de iesire pentru scriere
    with open(output_file, "w") as out_file:
        # Scriem header-ul în fisierul de iesire
        out_file.writelines(header)

        # Initializam dictionare pentru a stoca secventele si informatiile gene si CDS #M.I.
        gene_info_dict = {}
        cds_info_dict = {}

        # Iteram prin toate caracteristicile (features) din fisierul GenBank
        for feature in record.features:
            if feature.type == "gene":
                # Extragem informatiile despre gene
                gene_info = extract_feature_info(feature, record)
                gene_id = gene_info["gene_id"]
                gene_info_dict[gene_id] = gene_info
            elif feature.type == "CDS":
                # Extragem informatiile despre CDS
                cds_info = extract_feature_info(feature, record)
                gene_id = cds_info["gene_id"]
                cds_info_dict[gene_id] = cds_info

        # Scriem informatiile despre gene si CDS-uri în fisierul de iesire
        for gene_id, gene_info in gene_info_dict.items():
            # Scriem informatiile despre gene
            out_file.write(f"Gene: {gene_info['location']}\n")
            out_file.write(f" /gene={gene_id}\n")
            if "locus_tag" in gene_info:
                out_file.write(f" /locus_tag=\"{gene_info['locus_tag']}\"\n")
            if gene_info["note"]:
                out_file.write(f" /note=\"{gene_info['note']}\"\n")
            out_file.write("\n")

            # Scriem informatiile despre CDS daca exista
            if gene_id in cds_info_dict:
                cds_info = cds_info_dict[gene_id]
                out_file.write(f"CDS: {cds_info['location']}\n")
                out_file.write(f" /gene={gene_id}\n")
                if cds_info["protein_id"]:
                    out_file.write(f" /protein_id=\"{cds_info['protein_id']}\"\n")
                if cds_info["note"]:
                    out_file.write(f" /note=\"{cds_info['note']}\"\n")
                out_file.write("\n")

            # Combinam secventa genei si secventa CDS-ului (daca exista) si o scriem #M.I.
            gene_seq = gene_info["sequence"]
            cds_seq = cds_info_dict[gene_id]["sequence"] if gene_id in cds_info_dict else ""
            full_sequence = str(gene_seq) + str(cds_seq)
            out_file.write(f" /sequence:\n{full_sequence}\n\n")

def extract_feature_info(feature, record):
    # Extragem informatiile necesare pentru gene/CDS
    location = feature.location
    qualifiers = feature.qualifiers
    sequence = feature.extract(record.seq)

    # Obtinem identificatorii genei si alte note din calificative
    gene_id = qualifiers.get("gene", ["N/A"])[0]
    note = qualifiers.get("note", [""])[0]
    protein_id = qualifiers.get("protein_id", [""])[0]

    # Cream un dictionar cu informatiile relevante
    info = {
        "location": location,
        "gene_id": gene_id,
        "note": note,
        "protein_id": protein_id,
        "sequence": sequence
    }

    if "locus_tag" in qualifiers:
        info["locus_tag"] = qualifiers["locus_tag"][0]
    return info

if __name__ == "__main__":
    # Specificam fisierul de intrare si iesire
    input_file = r"C:\Users\tynu2\Desktop\pyt\Homo_sapiens.GRCh38.111.chromosome.19.dat"
    output_file = r"C:\Users\tynu2\Desktop\pyt\output.txt"
    # Apelam functia principala pentru extragerea genelor si CDS-urilor #M.I.
    extract_genes_and_cds(input_file, output_file)