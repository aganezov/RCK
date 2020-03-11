from setuptools import setup

import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import rck

setup(
    name="RCK",
    version=rck.version,
    author="Sergey Aganezov",
    author_email="aganezov@cs.jhu.edu",
    description="A tool for (R)econstruction of (C)ancer (K)aryotypes (both clone- and haplotype-specific)",
    license="MIT",
    keywords="RCK, rck, cancer, cancer genomics, cancer karyotypes, clonality, subclonality, copy number aberrations, breakpoints, structural variations, novel adjacencies",
    url="https://github.com/aganezov/rck",
    packages=["", "rck", "rck.core", "rck.utils", "rck.utils.scn", "rck.utils.adj"],
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "rck = rck.rck_run:main",
            "rck-scnt-x2rck = rck.utils.scn.rck_scnt_x2rck:main",
            "rck-scnt-process = rck.utils.scn.rck_scnt_process:main",
            "rck-scnt-rck2x = rck.utils.scn.rck_scnt_rck2x:main",
            "rck-scnt-stats = rck.utils.scn.rck_scnt_stats:main",
            "rck-scnb = rck.utils.scn.rck_scnb:main",
            "rck-adj-x2rck = rck.utils.adj.rck_adj_x2rck:main",
            "rck-adj-rck2x = rck.utils.adj.rck_adj_rck2x:main",
            "rck-adj-process = rck.utils.adj.rck_adj_process:main",
            "rck-adj-stats = rck.utils.adj.rck_adj_stats:main",
            "rck-adg-infer = rck.utils.adj.rck_adg_infer:main",
            "rck-adg-process = rck.utils.adj.rck_adg_process:main",
            "rck-adg-stats = rck.utils.adj.rck_adg_stats:main",
            "rck-input-refine = rck.utils.rck_input_refine:main",
            "rck-kar-graph = rck.utils.karyotype.rck_kar_graph:main",
            "rck-kar-stats = rck.utils.karyotype.rck_kar_stats:main",
        ]
    },
    install_requires=[
        "networkx>=2",
        "scipy",
        "pyvcf",
        "pysam",
        "sortedcontainers",
        "pandas",
        "gffutils",
    ]
)