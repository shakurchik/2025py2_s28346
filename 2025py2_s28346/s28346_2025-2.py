#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt

class NCBIRetriever:
    def __init__(self, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key
        self.webenv = None
        self.query_key = None

    def search_taxid(self, taxid):
        """Search for nucleotide records."""
        handle = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]", usehistory="y")
        result = Entrez.read(handle)
        self.webenv = result["WebEnv"]
        self.query_key = result["QueryKey"]
        return int(result["Count"])

    def fetch_records(self, max_records=1000):
        """Fetch GenBank records."""
        handle = Entrez.efetch(
            db="nucleotide",
            rettype="gb",
            retmode="text",
            retstart=0,
            retmax=max_records,
            webenv=self.webenv,
            query_key=self.query_key
        )
        return list(SeqIO.parse(handle, "genbank"))

    def filter_records(self, records, min_len, max_len):
        """Return only records within length range."""
        return [rec for rec in records if min_len <= len(rec.seq) <= max_len]

    def generate_csv(self, records, filename="records.csv"):
        """Generate CSV."""
        data = [{
            "Accession": rec.id,
            "Length": len(rec.seq),
            "Description": rec.description
        } for rec in records]
        df = pd.DataFrame(data)
        df.to_csv(filename, index=False)
        print(f"CSV saved to {filename}")
        return df

    def plot_lengths(self, df, filename="length_chart.png"):
        """Generate and save chart."""
        df_sorted = df.sort_values(by="Length", ascending=False)
        plt.figure(figsize=(12, 6))
        plt.plot(df_sorted["Accession"], df_sorted["Length"], marker='o')
        plt.xticks(rotation=90)
        plt.xlabel("Accession")
        plt.ylabel("Sequence Length")
        plt.title("GenBank Sequences by Length")
        plt.tight_layout()
        plt.savefig(filename)
        print(f"Chart saved to {filename}")


def main():

    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")
    taxid = input("Enter taxonomic ID (e.g. 9606 for human): ")
    min_len = int(input("Enter minimum sequence length: "))
    max_len = int(input("Enter maximum sequence length: "))

    retriever = NCBIRetriever(email, api_key)
    total = retriever.search_taxid(taxid)
    print(f"Found {total} records. Downloading and filtering...")

    all_records = retriever.fetch_records(max_records=min(1000, total))
    filtered = retriever.filter_records(all_records, min_len, max_len)
    print(f"{len(filtered)} records matched the length filter.")

    df = retriever.generate_csv(filtered)
    retriever.plot_lengths(df)


if __name__ == "__main__":
    main()
