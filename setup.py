from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="mordred-descriptor-calculator",
    version="0.1.0",
    author="377H-Miru",
    description="Professional Mordred and Conjugation Molecular Descriptor Calculator CLI",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/377H-Miru/mordred_descriptor_calculator",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=[
        "rdkit",
        "mordred",
        "networkx",
        "pandas",
        "numpy<2.0.0",
        "tqdm",
    ],
    extras_require={
        "dev": ["pytest"],
    },
    entry_points={
        "console_scripts": [
            "mordred-desc=mordred_descriptor_calculator.cli:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
)
