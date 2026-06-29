#!/usr/bin/env python
"""
Backward-compatibility wrapper for calculate_descriptors.py.
Delegates to mordred_descriptor_calculator.cli.main().
"""
from mordred_descriptor_calculator.cli import main

if __name__ == "__main__":
    main()
