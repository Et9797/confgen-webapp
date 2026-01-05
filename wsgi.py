#!/usr/bin/env python
import sys
import os

# Add app directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import Flask app
from app import app as application
