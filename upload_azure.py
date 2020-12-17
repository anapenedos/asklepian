import os
import sys
import argparse

from azure.storage.blob import BlobServiceClient
from azure.core.exceptions import ResourceNotFoundError, ResourceExistsError

azure_endpoint = os.environ.get('AZURE_END')
azure_saskey = os.environ.get('AZURE_SAS')

if not azure_endpoint or not azure_saskey:
    sys.stderr.write("[FAIL] AZURE_END or AZURE_SAS environment variable is unset.\n")
    sys.exit(1)

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--container", help="Base container", required=True)
parser.add_argument("-f", "--file", help="Local file (to upload or download)")
parser.add_argument("-b", "--blob", help="Remote file (to download)")
args = parser.parse_args()

# Connect to store and acquire container client
blob_service_client = BlobServiceClient(account_url=azure_endpoint, credential=azure_saskey)
container_name = args.container
try:
    container_client = blob_service_client.get_container_client(container_name)
except ResourceNotFoundError:
    container_client = blob_service_client.create_container(container_name, metadata=None, public_access=None)
    print("[NOTE] Container %s not found. Created it." % container_name)


if args.file and not args.blob:
    # If only a local file name provided, upload the blob
    if not os.path.isfile(args.file):
        sys.stderr.write("[FAIL] Cannot open %s" % args.file)
        sys.exit(2)

    remote_filename = os.path.basename(args.file)
    blob_client = blob_service_client.get_blob_client(container=container_name, blob=remote_filename)
    with open(args.file, "rb") as data:
        try:
            blob_client.upload_blob(data)
        except ResourceExistsError:
            sys.stderr.write("[FAIL] Remote blob named %s already exists. Refusing to overwrite.\n" % remote_filename)
            sys.exit(3)

elif args.blob:
    if args.file:
        download_to = args.file
    else:
        download_to = os.path.abspath(args.blob)

    blob_client = blob_service_client.get_blob_client(container=container_name, blob=args.blob)
    with open(download_to, "wb") as download_fh:
        try:
            download_fh.write(blob_client.download_blob().readall())
        except:
            sys.stderr.write("[FAIL] Error encountered downloading remote blob %s.\n" % args.blob)
            sys.exit(3)

# List all dirs and blobs
sys.stderr.write("[NOTE] Walking container %s:\n" % container_name)
blobs = container_client.walk_blobs()
for blob in blobs:
    sys.stderr.write("\t%s (%d)\n" % (blob.name, blob.size))
    #for k, v in blob.metadata.items():
    #    sys.stderr.write("\t\t%s\t\n" % (k, str(v)))

