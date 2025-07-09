from dataclasses import dataclass, field
import dataclasses
import os
from typing import Literal
import pysam
from pathlib import Path
from rnacloud_genome_reference.gtf import Feature, Exon, Intron, GTFHandler
from rnacloud_genome_reference.config import Config

import logging

logger = logging.getLogger(__name__)

@dataclass
class FeatureComparisonResult:
    primary_contig_transcript: str | None
    primary_contig_n_exons: int
    primary_contig_n_introns: int
    fix_contig_transcript: str | None
    fix_contig_n_exons: int
    fix_contig_n_introns: int
    n_exons_equal: bool
    n_introns_equal: bool
    sequences_unequal_n_exons: int
    sequences_unequal_n_introns: int
    splice_sites_unequal_n: int
    primary_exon_lengths: str | None = field(default=None, repr=False)
    primary_intron_lengths: str | None = field(default=None, repr=False)
    fix_exon_lengths: str | None = field(default=None, repr=False)
    fix_intron_lengths: str | None = field(default=None, repr=False)
    discordant_exon_numbering: bool | None = field(default=None, repr=False)
    comparison_status: Literal["Identical",
                               "Different - No. of exons or introns differ",
                               "Different - Sequences differ",
                               "Different - Unknown reason",
                               "Different - Splice-site sequences differ",
                               "Different - Exon numbering is discordant",
                               "Not comparable",
                               "Unknown"] = field(default="Unknown", init=False)

    def __post_init__(self): 
        if (self.n_exons_equal and self.n_introns_equal and self.sequences_unequal_n_exons == 0 and self.sequences_unequal_n_introns == 0):
            self.comparison_status = "Identical"
        
        if not self.n_exons_equal or not self.n_introns_equal:
            self.comparison_status = "Different - No. of exons or introns differ"
        
        if self.n_exons_equal and self.n_introns_equal and (self.sequences_unequal_n_exons > 0 or self.sequences_unequal_n_introns > 0):
            self.comparison_status = "Different - Sequences differ"
        
        if self.splice_sites_unequal_n > 0:
            self.comparison_status = "Different - Splice-site sequences differ"

        if self.fix_contig_transcript is not None and self.primary_contig_transcript is not None and self.discordant_exon_numbering is True:
            self.comparison_status = "Different - Exon numbering is discordant"

        if self.primary_contig_transcript is not None and self.fix_contig_transcript is None:
            self.comparison_status = "Not comparable"
        
class FeatureSequenceHelper:
    def __init__(self, fasta: str):
        self.fasta = pysam.FastaFile(fasta)

    def get_seq_for_feature(self, features: list[Feature]) -> list[Feature]:
        for feature in features:
            if feature.sequence is None:
                try:
                    feature.sequence = self.fasta.fetch(
                        feature.chromosome, 
                        feature.start - 1,  # pysam is 0-based
                        feature.end
                    ).upper()
                except ValueError as e:
                    logger.error(f"Error fetching sequence for {feature.chromosome}:{feature.start}-{feature.end}: {e}")
                    raise
        return features
    
    @staticmethod
    def get_feature_seq_lengths(features: list[Feature]) -> str | None:
        """
        Get a string representation of the lengths of the features.
        :param features: List of features to get lengths for.
        :return: String representation of feature lengths.
        """
        if len(features) == 0:
            logger.warning("No features provided. Cannot compute lengths.")
            return None
        
        if features[0].sequence is None:
            logger.error("Features do not have sequences. Cannot compute lengths.")
            raise ValueError("Features do not have sequences.")
        else:
            return ";".join(str(len(feature.sequence)) for feature in features) # type: ignore

class FeatureComparator:
    def __init__(self, gtf_file_path: str, fasta_file_path: str):
        self.gtf_file_path = gtf_file_path
        self.fasta_file_path = fasta_file_path

    @staticmethod
    def flag_discordant_exon_numbering(primary_features: list[Exon], fix_features: list[Exon]) -> bool | None:
        if len(primary_features) != len(fix_features):
            return None
        
        for primary_feature, fix_feature in zip(primary_features, fix_features):
            logger.debug("🤔 Comparing exons:")
            logger.debug(f"  Primary Feature: {primary_feature}")
            logger.debug(f"  Fix Feature:     {fix_feature}")
            if primary_feature.exon_no != fix_feature.exon_no:
                logger.debug("    ‼️ Exon numbers are discordant.")
                return True
            
        return False

    @staticmethod
    def compare_sequences(primary_features: list[Feature], fix_features: list[Feature]) -> int:
        """
        Compare sequences of primary and fix features.
        Returns the number of features with unequal sequences.
        If the lengths of the lists are not equal, returns -1.
        :param primary_features: List of primary features to compare.
        :param fix_features: List of fix features to compare.
        :return: Number of features with unequal sequences or -1 if lengths are not equal
        """
        if len(primary_features) != len(fix_features):
            return -1
        
        features_unequal = 0

        for primary_feature, fix_feature in zip(primary_features, fix_features):
            logger.debug("🤔 Comparing sequences:")
            logger.debug(f"  Primary Feature: {primary_feature}")
            logger.debug(f"  Fix Feature:     {fix_feature}")
            if primary_feature.sequence != fix_feature.sequence:
                logger.debug("    ‼️ Sequences are unequal.")
                features_unequal += 1

        return features_unequal
    
    @staticmethod
    def compare_splice_site_motifs(primary_features: list[Intron], fix_features: list[Intron]) -> int:
        splice_sites_unequal = 0

        if len(primary_features) != len(fix_features):
            return -1
        
        for primary_feature, fix_feature in zip(primary_features, fix_features):
            logger.debug("🤔 Comparing splice site motifs:")
            logger.debug(f"  Primary Feature: {primary_feature}")
            logger.debug(f"  Fix Feature:     {fix_feature}")
            if primary_feature.sequence is not None and fix_feature.sequence is not None:
                if primary_feature.sequence[0:2] != fix_feature.sequence[0:2]:
                    logger.debug(f"    ‼️ Splice site motif is unequal. Primary: {primary_feature.sequence[0:2]}, Fix: {fix_feature.sequence[0:2]}")
                    splice_sites_unequal += 1

                if primary_feature.sequence[-2:] != fix_feature.sequence[-2:]:
                    logger.debug(f"    ‼️ Splice site motif is unequal. Primary: {primary_feature.sequence[-2:]}, Fix: {fix_feature.sequence[-2:]}")
                    splice_sites_unequal += 1
        
        return splice_sites_unequal

    def compare_features(self, 
                         primary_chromosome: str, 
                         primary_start: int, 
                         primary_end: int,
                         fix_chromosome: str, 
                         fix_start: int, 
                         fix_end: int, 
                         entrez_gene_id: int):  
        gtf_handler = GTFHandler(self.gtf_file_path)
        logger.info(f"Comparing regions: Primary {primary_chromosome}:{primary_start}-{primary_end} with Fix {fix_chromosome}:{fix_start}-{fix_end} for Gene ID: {entrez_gene_id}")

        primary_transcript = self.obtain_transcript(primary_chromosome, primary_start, primary_end, entrez_gene_id, gtf_handler)

        primary_exons = gtf_handler.get_exons_by_transcript(
            chromosome=primary_chromosome,
            start=primary_start,
            end=primary_end,
            transcript_id=primary_transcript # type: ignore
        )

        primary_introns = gtf_handler.derive_introns_from_exons(exons=primary_exons)

        fix_exons = gtf_handler.get_exons_by_transcript(
            chromosome=fix_chromosome,
            start=fix_start,
            end=fix_end,
            transcript_id=primary_transcript # type: ignore
        )

        # If no exons are found in the fix region, we log a warning and set the fix transcript to None.
        # This is to ensure that we do not proceed with an empty list of exons,
        # which would lead to errors in subsequent processing.
        # This can occur if the primary transcript does not have a corresponding fix region in the GTF file.
        if len(fix_exons) == 0:
            logger.warning(f"No exons found in fix region {fix_chromosome}:{fix_start}-{fix_end} for transcript {primary_transcript}.")
            fix_transcript = None
            fix_introns = []
        else:
            fix_transcript = primary_transcript
            fix_introns = gtf_handler.derive_introns_from_exons(exons=fix_exons)

        n_exons_equal = len(primary_exons) == len(fix_exons)
        n_introns_equal = len(primary_introns) == len(fix_introns)

        logger.debug(f"Exons equal: {n_exons_equal}, Introns equal: {n_introns_equal}")

        feature_sequence_helper = FeatureSequenceHelper(self.fasta_file_path)
        primary_exons = feature_sequence_helper.get_seq_for_feature(primary_exons)
        primary_introns = feature_sequence_helper.get_seq_for_feature(primary_introns)

        fix_exons = feature_sequence_helper.get_seq_for_feature(fix_exons)
        fix_introns = feature_sequence_helper.get_seq_for_feature(fix_introns)

        primary_exons_lengths = FeatureSequenceHelper.get_feature_seq_lengths(primary_exons)
        primary_introns_lengths = FeatureSequenceHelper.get_feature_seq_lengths(primary_introns)
        fix_exons_lengths = FeatureSequenceHelper.get_feature_seq_lengths(fix_exons)
        fix_introns_lengths = FeatureSequenceHelper.get_feature_seq_lengths(fix_introns)

        sequences_unequal_n_exons = FeatureComparator.compare_sequences(primary_exons, fix_exons)
        sequences_unequal_n_introns = FeatureComparator.compare_sequences(primary_introns, fix_introns)
        splice_sites_unequal_n = FeatureComparator.compare_splice_site_motifs(primary_introns, fix_introns)

        discordant_exon_numbering = FeatureComparator.flag_discordant_exon_numbering(primary_exons, fix_exons)

        logger.debug(f"No. of exons with unequal sequences: {sequences_unequal_n_exons}")
        logger.debug(f"No. of introns with unequal sequences: {sequences_unequal_n_introns}")

        result = FeatureComparisonResult(
            primary_contig_transcript=primary_transcript,
            primary_contig_n_exons=len(primary_exons),
            primary_contig_n_introns=len(primary_introns),
            fix_contig_transcript=fix_transcript,
            fix_contig_n_exons=len(fix_exons),
            fix_contig_n_introns=len(fix_introns),
            n_exons_equal=n_exons_equal,
            n_introns_equal=n_introns_equal,
            sequences_unequal_n_exons=sequences_unequal_n_exons,
            sequences_unequal_n_introns=sequences_unequal_n_introns,
            splice_sites_unequal_n=splice_sites_unequal_n,
            primary_exon_lengths=primary_exons_lengths,
            primary_intron_lengths=primary_introns_lengths,
            fix_exon_lengths=fix_exons_lengths,
            fix_intron_lengths=fix_introns_lengths,
            discordant_exon_numbering=discordant_exon_numbering
        )

        return dataclasses.asdict(result)

    def obtain_transcript(self, primary_chromosome, primary_start, primary_end, entrez_gene_id, gtf_handler):
        logger.info(f"Fetching transcript for Entrez Gene ID: {entrez_gene_id} in primary region {primary_chromosome}:{primary_start}-{primary_end}")

        transcript = gtf_handler.get_transcript_for_gene(
            primary_chromosome, 
            primary_start, 
            primary_end, 
            entrez_gene_id,
            True
        )

        if transcript is None:
            logger.warning(f"No MANE transcript found for Entrez Gene ID: {entrez_gene_id} in primary region {primary_chromosome}:{primary_start}-{primary_end}. Attempting to fetch non-MANE transcript.")
            transcript = gtf_handler.get_transcript_for_gene(
                primary_chromosome, 
                primary_start, 
                primary_end, 
                entrez_gene_id,
                False
            )

        if transcript is None:
            logger.warning(f"Could not find any transcript for Entrez Gene ID: {entrez_gene_id} in region {primary_chromosome}:{primary_start}-{primary_end}.")
            
        return transcript
    