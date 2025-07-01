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
    comparison_status: Literal["Identical",
                               "Different - No. of exons or introns differ",
                               "Different - Sequences differ",
                                "Different - Unknown reason",
                               "Not comparable"] = field(init=False)

    def __post_init__(self):
        if self.fix_contig_transcript is None or self.primary_contig_transcript is None:
            self.comparison_status = "Not comparable"
        elif (self.n_exons_equal and self.n_introns_equal and self.sequences_unequal_n_exons == 0 and self.sequences_unequal_n_introns == 0):
            self.comparison_status = "Identical"
        elif not self.n_exons_equal or not self.n_introns_equal:
            self.comparison_status = "Different - No. of exons or introns differ"
        elif self.n_exons_equal and self.n_introns_equal and (self.sequences_unequal_n_exons > 0 or self.sequences_unequal_n_introns > 0):
            self.comparison_status = "Different - Sequences differ"
        else:
            self.comparison_status = "Different - Unknown reason"
        
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

class FeatureComparator:
    def __init__(self, gtf_file_path: str, fasta_file_path: str):
        self.gtf_file_path = gtf_file_path
        self.fasta_file_path = fasta_file_path

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
        fix_exons = feature_sequence_helper.get_seq_for_feature(fix_exons)

        primary_introns = feature_sequence_helper.get_seq_for_feature(primary_introns)
        fix_introns = feature_sequence_helper.get_seq_for_feature(fix_introns)

        sequences_unequal_n_exons = FeatureComparator.compare_sequences(primary_exons, fix_exons)
        sequences_unequal_n_introns = FeatureComparator.compare_sequences(primary_introns, fix_introns)

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
            sequences_unequal_n_introns=sequences_unequal_n_introns
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
    