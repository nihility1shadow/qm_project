"""Compare compact reference reachability/order to the original graph oracle."""
from pathlib import Path
import argparse, subprocess, tempfile
root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--cxx", default="g++")
args = parser.parse_args()
source = (root/"ahm-mb-sep.cpp").read_text()
start = source.index("double sepmb_binom(")
stop = source.index("double sepmb_kondo_degeneracy(",start)
prefix = '#include <algorithm>\n#include <cmath>\n#include <iostream>\n#include <vector>\n#include <set>\n#include <unordered_map>\n#include <unordered_set>\nusing namespace std;\n'
main = 'int main() {\n  int cases=0;\n  for(int norb=3;norb<=12;norb++) for(int nel=1;nel<norb;nel++)\n    for(int mode=0;mode<2;mode++) for(int distance=0;distance<=4;distance++) {\n      set<int> initial;\n      for(int i=0;i<nel;i++) initial.insert(mode ? norb-1-i : i);\n      const auto original=sepmb_generate_reference_states(norb,initial,distance);\n      const auto compact=sepmb_generate_reference_masks(norb,initial,distance);\n      if(original.size()!=compact.size()) return 1;\n      if(original.size()!=sepmb_reference_state_count(norb,nel,initial.count(0),distance)) return 2;\n      for(size_t i=0;i<compact.size();i++)\n        if(sepmb_state_mask(original[i])!=compact[i]) return 3;\n      cases++;\n    }\n  for(const set<int> initial : {set<int>{0,63},set<int>{62,63}}) {\n    const auto original=sepmb_generate_reference_states(64,initial,0);\n    const auto compact=sepmb_generate_reference_masks(64,initial,0);\n    if(original.size()!=compact.size()) return 4;\n    for(size_t i=0;i<compact.size();i++)\n      if(sepmb_state_mask(original[i])!=compact[i]) return 5;\n    cases++;\n  }\n  cout << cases << " reference graphs and ordering checks passed\\n";\n}\n'
with tempfile.TemporaryDirectory(prefix="qm-reference-") as directory:
    directory = Path(directory)
    cpp = directory/"reference_test.cpp"
    exe = directory/"reference_test.exe"
    cpp.write_text(prefix+source[start:stop]+main)
    subprocess.run([args.cxx,"-O2","-std=c++17",str(cpp),"-o",str(exe)],check=True,timeout=120)
    subprocess.run([str(exe)],check=True,timeout=120)
