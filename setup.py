import setuptools

setuptools.setup(
    name='treeflow-paper',
    version='0.0.1',
    packages=setuptools.find_packages(),
    install_requires=[
        'treeflow',
        'pyyaml',
        'biopython',
        'jinja2',
        'dendropy',
        'pandas'
    ]
)
