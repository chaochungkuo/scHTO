from setuptools import setup, find_packages
from Cython.Build import cythonize

setup(
    name="schto",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "PyYAML",
        "pytest",
        "modin[all]",
        "numpy",
        "ray"
    ],
    entry_points={
        "console_scripts": [
            "schto=src.main:main"
        ]
    },
    ext_modules = cythonize("src/fastq_chunk_processor.pyx",
                            compiler_directives={'language_level' : "3"})
)