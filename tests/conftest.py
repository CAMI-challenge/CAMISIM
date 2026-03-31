"""Shared fixtures and path setup for CAMISIM tests."""

import sys
import os

# Add script directories to Python path so we can import modules
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "pipelines", "shared", "scripts"))
sys.path.insert(0, os.path.join(ROOT, "pipelines", "metagenomic", "scripts"))
