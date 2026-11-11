import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent))  # import recipe
sys.path.insert(0, str(HERE))         # fake_scorer (also importable from spawned workers)
