# add pyPANAIR to python path
from pathlib import Path
import sys
sys.path.append(str(Path.cwd().parent.parent))

from postprocess import ffmf_converter

if __name__ == '__main__':
    ffmf_converter.write_ffmf()
