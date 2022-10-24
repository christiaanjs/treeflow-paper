import setuptools

setuptools.setup(
    name="treeflow-paper",
    version="0.0.1",
    packages=setuptools.find_packages(),
    install_requires=[
        "treeflow",
        "pyyaml",
        "jinja2",
        "dendropy",
        "pandas",
        "click",
        "click-config-file",
        "importlib",
        "snakemake",
        "papermill",
        "biopython",
        "matplotlib",
        "humanfriendly",
    ],
    entry_points="""
        [console_scripts]
        treeflow_pipeline=treeflow_pipeline.cli:cli
    """,
    package_data={"": ["*.smk"]},
)
