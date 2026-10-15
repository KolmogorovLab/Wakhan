#!/usr/bin/env python3

import argparse
import logging
import zipfile
from pathlib import Path

logger = logging.getLogger()

def create_hiscanner_plot_data_archives(
    out_dir,
    genome_name,
    breakpoints,
    centromere_bed,
):
    """
    Create one HiScanner plots-data ZIP archive for every Wakhan solution.

    Expected directory structure:

        out_dir/
        ├── solution_1 -> solution_<...>
        ├── solution_2 -> solution_<...>
        ├── solution_rank_1 -> solution_<...>
        ├── solution_rank_2 -> solution_<...>
        ├── solutions_ranks.tsv
        │
        ├── coverage_data/
        │   ├── phase_corrected_coverage.csv
        │   ├── baf.csv
        │   └── <genome-name>_loh_segments.csv   # optional
        │
        ├── solution_<...>/
        │   ├── integer_profile.bed
        │   ├── subclonal_profile.bed
        │   └── genes/                          # optional
        │       └── genes_copynumber_states.bed
        └── ...

    Creates (inside each solution's own directory):

        solution_<...>/HiScanner_plots_data.zip
        ...

    Each ZIP contains:

        integer_profile.bed
        subclonal_profile.bed
        genes_copynumber_states.bed      # if available
        phase_corrected_coverage.csv
        baf.csv
        <genome-name>_loh_segments.csv   # if available
        <centromere-bed-name>
        severus_somatic.vcf              # if --breakpoints was provided
        solutions_ranks.tsv
    """

    out_dir = Path(out_dir).resolve()
    breakpoints = Path(breakpoints).resolve() if breakpoints else None
    centromere_bed = Path(centromere_bed).resolve()

    genome_name = str(genome_name)

    # ------------------------------------------------------------------
    # Validate output directory
    # ------------------------------------------------------------------

    if not out_dir.exists():
        raise FileNotFoundError(
            f"Working/output directory does not exist: {out_dir}"
        )

    if not out_dir.is_dir():
        raise NotADirectoryError(
            f"Output path is not a directory: {out_dir}"
        )

    # ------------------------------------------------------------------
    # Validate common files
    # ------------------------------------------------------------------

    if breakpoints is not None and not breakpoints.exists():
        raise FileNotFoundError(
            f"Breakpoints file does not exist: {breakpoints}"
        )

    if breakpoints is None:
        logger.info(
            "No --breakpoints provided: severus_somatic.vcf will "
            "be omitted from the archive(s)"
        )

    if not centromere_bed.exists():
        raise FileNotFoundError(
            f"Centromere BED file does not exist: {centromere_bed}"
        )

    # ------------------------------------------------------------------
    # solutions_ranks.tsv
    # ------------------------------------------------------------------

    solutions_ranks = out_dir / "solutions_ranks.tsv"

    if not solutions_ranks.exists():
        raise FileNotFoundError(
            f"Missing solutions ranks file: {solutions_ranks}"
        )

    # ------------------------------------------------------------------
    # Coverage files
    # ------------------------------------------------------------------

    coverage_dir = out_dir / "coverage_data"

    if not coverage_dir.exists():
        raise FileNotFoundError(
            f"Coverage directory does not exist: {coverage_dir}"
        )

    phase_corrected_coverage = (
        coverage_dir / "phase_corrected_coverage.csv"
    )

    baf = coverage_dir / "baf.csv"

    if not phase_corrected_coverage.exists():
        raise FileNotFoundError(
            f"Missing required coverage file: "
            f"{phase_corrected_coverage}"
        )

    if not baf.exists():
        raise FileNotFoundError(
            f"Missing required BAF file: {baf}"
        )

    # ------------------------------------------------------------------
    # Optional LOH segments file
    # ------------------------------------------------------------------

    loh_segments = (
        coverage_dir / f"{genome_name}_loh_segments.csv"
    )

    if loh_segments.exists():
        logger.info(
            "Found optional LOH segments file: %s",
            loh_segments,
        )
    else:
        logger.info(
            "Optional LOH segments file not found: %s",
            loh_segments,
        )

    # ------------------------------------------------------------------
    # Find solution symlinks
    #
    # Supported:
    #
    #   solution_1
    #   solution_2
    #
    # and:
    #
    #   solution_rank_1
    #   solution_rank_2
    # ------------------------------------------------------------------

    solution_links = {}

    for path in out_dir.iterdir():

        if not path.is_symlink():
            continue

        name = path.name

        # solution_1, solution_2, ...
        if name.startswith("solution_"):
            suffix = name[len("solution_"):]

            if suffix.isdigit():
                rank = int(suffix)
                solution_links[rank] = path
                continue

        # solution_rank_1, solution_rank_2, ...
        if name.startswith("solution_rank_"):
            suffix = name[len("solution_rank_"):]

            if suffix.isdigit():
                rank = int(suffix)
                solution_links[rank] = path

    if not solution_links:
        logger.warning(
            "No solution symlinks found in %s",
            out_dir,
        )
        return []

    solutions = sorted(
        solution_links.items(),
        key=lambda item: item[0],
    )

    logger.info(
        "Found %d solution(s)",
        len(solutions),
    )

    # ------------------------------------------------------------------
    # Create ZIP for every solution
    # ------------------------------------------------------------------

    created_archives = []

    for solution_number, solution_link in solutions:

        solution_dir = solution_link.resolve()

        logger.info(
            "Processing solution %d: %s -> %s",
            solution_number,
            solution_link.name,
            solution_dir,
        )

        if not solution_dir.exists():
            raise FileNotFoundError(
                f"Solution {solution_number} points to a "
                f"non-existent directory: {solution_dir}"
            )

        if not solution_dir.is_dir():
            raise NotADirectoryError(
                f"Solution {solution_number} is not a directory: "
                f"{solution_dir}"
            )

        # --------------------------------------------------------------
        # Solution-specific files
        # --------------------------------------------------------------

        integer_profile = (
            solution_dir / "integer_profile.bed"
        )

        subclonal_profile = (
            solution_dir / "subclonal_profile.bed"
        )

        # Optional: genes/genes_copynumber_states.bed
        #
        # Not produced when haplotype phase correction is
        # bypassed (--phaseblock-flipping-disable, implied by
        # --without-phasing), since gene copy-number states are
        # only meaningful for phase-corrected segments.
        genes_copynumber_states = (
            solution_dir
            / "genes"
            / "genes_copynumber_states.bed"
        )

        if not integer_profile.exists():
            raise FileNotFoundError(
                f"Missing integer profile for solution "
                f"{solution_number}: {integer_profile}"
            )

        if not subclonal_profile.exists():
            raise FileNotFoundError(
                f"Missing subclonal profile for solution "
                f"{solution_number}: {subclonal_profile}"
            )

        if genes_copynumber_states.exists():
            logger.info(
                "Found optional gene copy-number states file: %s",
                genes_copynumber_states,
            )
        else:
            logger.info(
                "Optional gene copy-number states file not found: %s",
                genes_copynumber_states,
            )

        # --------------------------------------------------------------
        # ZIP filename
        # --------------------------------------------------------------

        zip_name = "HiScanner_plots_data.zip"

        zip_path = solution_dir / zip_name

        # Remove existing ZIP.
        if zip_path.exists():
            logger.info(
                "Removing existing archive: %s",
                zip_path,
            )
            zip_path.unlink()

        logger.info(
            "Creating: %s",
            zip_path,
        )

        # --------------------------------------------------------------
        # Create archive
        # --------------------------------------------------------------

        with zipfile.ZipFile(
            zip_path,
            mode="w",
            compression=zipfile.ZIP_DEFLATED,
            compresslevel=6,
        ) as zf:

            # ----------------------------------------------------------
            # Solution-specific files
            # ----------------------------------------------------------

            zf.write(
                integer_profile,
                arcname="integer_profile.bed",
            )

            zf.write(
                subclonal_profile,
                arcname="subclonal_profile.bed",
            )

            if genes_copynumber_states.exists():
                zf.write(
                    genes_copynumber_states,
                    arcname="genes_copynumber_states.bed",
                )

                logger.info(
                    "Added gene copy-number states: %s",
                    genes_copynumber_states.name,
                )

            # ----------------------------------------------------------
            # Coverage files
            # ----------------------------------------------------------

            zf.write(
                phase_corrected_coverage,
                arcname="phase_corrected_coverage.csv",
            )

            zf.write(
                baf,
                arcname="baf.csv",
            )

            # ----------------------------------------------------------
            # Optional LOH segments
            # ----------------------------------------------------------

            if loh_segments.exists():
                zf.write(
                    loh_segments,
                    arcname=loh_segments.name,
                )

                logger.info(
                    "Added LOH segments: %s",
                    loh_segments.name,
                )

            # ----------------------------------------------------------
            # Centromere BED
            # ----------------------------------------------------------

            zf.write(
                centromere_bed,
                arcname=centromere_bed.name,
            )

            # ----------------------------------------------------------
            # SeveRus somatic VCF
            #
            # Regardless of the original --breakpoints filename,
            # store it as severus_somatic.vcf.
            # ----------------------------------------------------------

            if breakpoints is not None:
                zf.write(
                    breakpoints,
                    arcname="severus_somatic.vcf",
                )

            # ----------------------------------------------------------
            # Wakhan solution ranking
            # ----------------------------------------------------------

            zf.write(
                solutions_ranks,
                arcname="solutions_ranks.tsv",
            )

        created_archives.append(zip_path)

        logger.info(
            "Created archive: %s",
            zip_path,
        )

    logger.info(
        "Finished creating %d HiScanner archive(s)",
        len(created_archives),
    )

    return created_archives
