#!/usr/bin/env python3

import boto3
import logging
import argparse

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s : %(message)s", force=True)

def download_file_from_s3(bucket, prefix, s3_client, output_path):

    try:
        s3_client.download_file(bucket, prefix, output_path)
        logging.info(f"Successfully downloaded s3://{bucket}/{prefix} to {output_path}")

    except Exception as e:
        logging.error(f"Error downloading file: {e}")

def parse_s3_uri(s3_uri):

    logging.debug("Get bucket and prefix")
    bucket_prefix = s3_uri.strip('s3://').strip('S3://').split('/')
    bucket = bucket_prefix[0]
    prefix = '/'.join(bucket_prefix[1:])

    logging.info(f"Parsed bucket: {bucket}, prefix: {prefix}")

    return bucket, prefix

def process_mapping_file(bucket, prefix, output_path):

    logging.debug("Initialize AWS Client")
    s3_client = boto3.client('s3')

    download_file_from_s3(bucket, prefix, s3_client, output_path)

    df = pd.read_csv(output_path, sep='\t')

    return df

def get_uri_from_map(df, accession):

    logging.debug("Find filename for given accession")

    try:
        row = df[df['accession'] == accession]
        if not row.empty:
            logging.info(f"Found uri: {row['filename'].iloc[0]} for accession: {accession}")
            return row['filename'].iloc[0]
        logging.warning(f"Accession {accession} not found in mapping file.")
        return None

    except Exception as e:
        logging.error(f"Error reading mapping file: {e}")
        return None

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog='Pull Accession from S3',
        description="Uses the gambit closest accession to pull the full S3 URI from the genome files."
    )
    parser.add_argument(
        "-m", "--mapping-file",
        required=True,
        help="S3 URI to TSV mapping file"
    )
    parser.add_argument(
        "-o", "--output-path",
        required=True,
        help="Local path to save the downloaded genome file"
    )
    parser.add_argument(
        "-a", "--accession",
        required=True,
        help="Accession number to find (e.g., GCF_12345)"
    )

    logging.debug("Run parser to call arguments downstream")
    args = parser.parse_args()

    bucket, prefix = parse_s3_uri(args.mapping_file)

    df = process_mapping_file(bucket, prefix, args.output_path)

    uri = get_uri_from_map(df, args.accession)

    uri_bucket, uri_prefix = parse_s3_uri(uri)

    download_file_from_s3(uri_bucket, uri_prefix, boto3.client('s3'), f"closest_accession_{args.accession}.fna")