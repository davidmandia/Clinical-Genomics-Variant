import sqlite3
import pysam
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

def get_transcripts_with_variants(database_path):
    # Connect to the SQLite3 database
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()
    
    # Query to get the unique transcript IDs with variants
    cursor.execute("SELECT DISTINCT transcript_ref FROM variants")
    transcripts_with_variants = {row[0] for row in cursor.fetchall()}
    
    conn.close()
    return transcripts_with_variants

def extract_transcript_lengths(sam_filename, valid_transcripts=None):
    lengths = []
    samfile = pysam.AlignmentFile(sam_filename, "r")

    for i, read in enumerate(samfile.fetch()):
        if valid_transcripts is None or read.query_name in valid_transcripts:  # Check if the read's transcript is in the list of valid transcripts
            aligned_length = 0
            if read.cigartuples is not None:
                # Sum the lengths of the 'M', '=', and 'X' operations
                aligned_length = sum([length for operation, length in read.cigartuples if operation in [0, 7, 8]])  # 0 = M, 7 = '=', 8 = 'X'
            
            lengths.append(aligned_length)
    
    samfile.close()
    return lengths

# Path to your SQLite3 database
database_path_38 = "output/database/GRCh38_indels_variant.db"
database_path_37 = "output/database/GRCh37_indels_variant.db"

# Get the list of transcripts with variants
transcripts_with_variants_38 = get_transcripts_with_variants(database_path_38)
transcripts_with_variants_37 = get_transcripts_with_variants(database_path_37)

# Process the SAM files for GRCh37 and GRCh38, filtering only those transcripts with variants
lengths_grch37_variants = extract_transcript_lengths("input/GRCh37.sam", transcripts_with_variants_37)
lengths_grch38_variants = extract_transcript_lengths("input/GRCh38.sam", transcripts_with_variants_38)

# Plotting the distributions using density and line plot
plt.figure(figsize=(10, 6))

# Plot density for GRCh37
counts_grch37, bins_grch37 = np.histogram(lengths_grch37_variants, bins=50, density=True)
plt.plot(bins_grch37[:-1], counts_grch37, label="GRCh37", linestyle='-', linewidth=2)

# Plot density for GRCh38
counts_grch38, bins_grch38 = np.histogram(lengths_grch38_variants, bins=50, density=True)
plt.plot(bins_grch38[:-1], counts_grch38, label="GRCh38", linestyle='-', linewidth=2)

# Add titles and labels
plt.title("Transcript Length Density Distribution for Transcripts with Variants")
plt.xlabel("Transcript Length")
plt.ylabel("Density")
plt.legend()

# Save the plot
plt.savefig("length_density_with_variants_analysis.png")

# Show the plot
plt.show()

# Perform t-test on transcripts with variants
t_stat_variants, p_value_variants = stats.ttest_ind(lengths_grch37_variants, lengths_grch38_variants, equal_var=False)

# Calculate means and standard deviations for transcripts with variants
mean_grch37_variants = np.mean(lengths_grch37_variants)
std_grch37_variants = np.std(lengths_grch37_variants)
mean_grch38_variants = np.mean(lengths_grch38_variants)
std_grch38_variants = np.std(lengths_grch38_variants)

# Print the results of the t-test for transcripts with variants
print("T-test for transcripts with variants between GRCh37 and GRCh38:")
print(f"T-statistic: {t_stat_variants:.4f}, P-value: {p_value_variants:.4e}")

# Print the means and standard deviations for transcripts with variants
print("\nStatistics for GRCh37 (with variants):")
print(f"Mean: {mean_grch37_variants:.2f}, Standard Deviation: {std_grch37_variants:.2f}")

print("\nStatistics for GRCh38 (with variants):")
print(f"Mean: {mean_grch38_variants:.2f}, Standard Deviation: {std_grch38_variants:.2f}")

# Now perform the same analysis for all transcripts (not just those with variants)
# Extract lengths for all transcripts
lengths_grch37_all = extract_transcript_lengths("input/GRCh37.sam")
lengths_grch38_all = extract_transcript_lengths("input/GRCh38.sam")

# Perform t-test on all transcripts
t_stat_all, p_value_all = stats.ttest_ind(lengths_grch37_all, lengths_grch38_all, equal_var=False)

# Calculate means and standard deviations for all transcripts
mean_grch37_all = np.mean(lengths_grch37_all)
std_grch37_all = np.std(lengths_grch37_all)
mean_grch38_all = np.mean(lengths_grch38_all)
std_grch38_all = np.std(lengths_grch38_all)

# Print the results of the t-test for all transcripts
print("\nT-test for all transcripts between GRCh37 and GRCh38:")
print(f"T-statistic: {t_stat_all:.4f}, P-value: {p_value_all:.4e}")

# Print the means and standard deviations for all transcripts
print("\nStatistics for GRCh37 (all transcripts):")
print(f"Mean: {mean_grch37_all:.2f}, Standard Deviation: {std_grch37_all:.2f}")

print("\nStatistics for GRCh38 (all transcripts):")
print(f"Mean: {mean_grch38_all:.2f}, Standard Deviation: {std_grch38_all:.2f}")
