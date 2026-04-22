from setuptools import setup, find_packages

setup(
    name="eskape_metquest",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "cobra",        # For ID mapping and model parsing
        "pandas",       # For result exporting
        "networkx",     # Likely used by your metquest_clone
    ],
    entry_points={
        'console_scripts': [
            'run-eskape=run_analysis:main',
            'run-metquest=run:main',
        ],
    },
)