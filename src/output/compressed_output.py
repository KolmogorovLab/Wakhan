#!/usr/bin/env python3

import argparse
import logging
import zipfile
from pathlib import Path


logger = logging.getLogger("wakhan_hiscanner_zip")


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
        ├── solution_rank_1 -> solution_<...>   # also supported
        ├── solution_rank_2 -> solution_<...>
        │
        ├── coverage_data/
        │   ├── phase_corrected_coverage.csv
        │   ├── baf.csv
        │   └── <genome-name>_loh_segments.csv   # optional
        │
        ├── solution_<...>/
        │   ├── integer_profile.bed
        │   └── subclonal_profile.bed
        └── ...

    Creates:

        <genome-name>_solution_1_HiScanner_plots_data.zip
        <genome-name>_solution_2_HiScanner_plots_data.zip
        ...

    Each ZIP contains:

        integer_profile.bed
        subclonal_profile.bed
        phase_corrected_coverage.csv
        baf.csv
        <genome-name>_loh_segments.csv   # if available
        <centromere-bed-name>
        severus_somatic.vcf
    """

    out_dir = Path(out_dir).resolve()
    breakpoints = Path(breakpoints).resolve()
    centromere_bed = Path(centromere_bed).resolve()

    genome_name = str(genome_name)

    if not out_dir.exists():
        raise FileNotFoundError(
            f"Working/output directory does not exist: {out_dir}"
        )

    # ------------------------------------------------------------------
    # Validate common input files
    # ------------------------------------------------------------------

    if not breakpoints.exists():
        raise FileNotFoundError(
            f"Breakpoints file does not exist: {breakpoints}"
        )

    if not centromere_bed.exists():
        raise FileNotFoundError(
            f"Centromere BED file does not exist: {centromere_bed}"
        )

    coverage_dir = out_dir / "coverage_data"

    if not coverage_dir.exists():
        raise FileNotFoundError(
            f"Coverage directory does not exist: {coverage_dir}"
        )

    phase_corrected_coverage = (
        coverage_dir / "phase_corrected_coverage.csv"
    )

    baf = coverage_dir / "baf.csv"

    # Optional file.
    #
    # Example:
    #   --genome-name H2009
    #
    # looks for:
    #   coverage_data/H2009_loh_segments.csv
    #
    loh_segments = (
        coverage_dir / f"{genome_name}_loh_segments.csv"
    )

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
    #
    # The solution symlink itself is not stored in the ZIP.
    # We resolve it and copy files from the actual target directory.
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
    # Report optional LOH file
    # ------------------------------------------------------------------

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

        # --------------------------------------------------------------
        # ZIP filename
        # --------------------------------------------------------------

        zip_name = (
            f"{genome_name}_solution_{solution_number}"
            f"_HiScanner_plots_data.zip"
        )

        zip_path = out_dir / zip_name

        # Remove existing ZIP so rerunning the command produces
        # a clean archive.
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

            # Solution-specific files
            zf.write(
                integer_profile,
                arcname="integer_profile.bed",
            )

            zf.write(
                subclonal_profile,
                arcname="subclonal_profile.bed",
            )

            # Common coverage files
            zf.write(
                phase_corrected_coverage,
                arcname="phase_corrected_coverage.csv",
            )

            zf.write(
                baf,
                arcname="baf.csv",
            )

            # ----------------------------------------------------------
            # Optional LOH segments file
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
            # store it in the archive as severus_somatic.vcf.
            # ----------------------------------------------------------

            zf.write(
                breakpoints,
                arcname="severus_somatic.vcf",
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
