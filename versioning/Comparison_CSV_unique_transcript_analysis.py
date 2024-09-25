import pandas as pd

# Load the two CSV files into DataFrames
file1 = 'versioning/results/highest_transcript_versions_GRCh37.csv'
file2 = 'versioning/results/highest_transcript_versions_GRCh38.csv'

df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Find transcript_base values present in one file but not in the other
bases_only_in_df1 = df1[~df1['transcript_base'].isin(df2['transcript_base'])]
bases_only_in_df2 = df2[~df2['transcript_base'].isin(df1['transcript_base'])]

# Find transcripts present in both files but with different versions
merged_df = pd.merge(df1, df2, on='transcript_base', suffixes=('_file1', '_file2'))
different_versions = merged_df[merged_df['transcript_version_file1'] != merged_df['transcript_version_file2']]

# Print the results for transcript_base differences
print("Transcript bases present only in the first file:")
print(bases_only_in_df1[['transcript_ref', 'transcript_base', 'transcript_version']])

print("\nTranscript bases present only in the second file:")
print(bases_only_in_df2[['transcript_ref', 'transcript_base', 'transcript_version']])

# Print the results for different versions
print("\nTranscripts with different versions between the two files:")
print(different_versions[['transcript_ref_file1', 'transcript_base', 'transcript_version_file1', 'transcript_ref_file2', 'transcript_version_file2']])

# Count the number of unique transcript_base values in each category
print(f"\nNumber of transcript bases only in the first file: {len(bases_only_in_df1)}")
print(f"Number of transcript bases only in the second file: {len(bases_only_in_df2)}")
print(f"Number of transcripts with different versions: {len(different_versions)}")
