"""Check cache round-trip and reject incomplete/corrupt/mismatched inputs."""
from pathlib import Path
import argparse,subprocess,tempfile
root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--cxx',default='g++');args=parser.parse_args()
source=(root/'ahm-mb-sep.cpp').read_text();start=source.index('unsigned long long sepmb_reference_checksum(');stop=source.index('struct SepmbReferenceTransitions {',start)
prefix='''#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cmath>
#include <iostream>
using namespace std;
'''
main=r'''
int main(int argc,char** argv) {
  if(argc!=2)return 1;
  double original[3][3]={{1,-0.0,1.e-310},{.1,.23456789012345678,4},{-.71,3,2}};
  double restored[3][3]={};double* rows[]={restored[0],restored[1],restored[2]};
  const string key="#SEP_MB_REFERENCE_KEY test";
  auto write=[&](int defect) {
    ofstream out(argv[1]);out<<key<<"\n# bounded reference\n";
    out<<scientific<<setprecision(16);
    unsigned long long sum=14695981039346656037ULL;
    for(int t=0;t<3;t++) {
      sum=sepmb_reference_checksum(sum,original[t],3);
      out<<(defect==3&&t==1?.3:t*.5);
      for(int c=0;c<3;c++) {
        if(defect==4&&t==1&&c==1)out<<" nan";
        else { double value=original[t][c]; if(defect==2&&t==1&&c==1)value+=.125; out<<' '<<value; }
      }
      out<<'\n';
    }
    if(defect!=1)out<<"#SEP_MB_REFERENCE_COMPLETE "<<sum<<'\n';
    if(defect==5)out<<"trailing data\n";
  };
  string error;write(0);
  if(!sepmb_read_reference_cache(argv[1],key,2,3,.5,rows,error))return 2;
  if(memcmp(original,restored,sizeof(original)))return 3;
  if(sepmb_read_reference_cache(argv[1],key+"wrong",2,3,.5,rows,error))return 4;
  for(int defect=1;defect<=5;defect++) {
    write(defect);
    if(sepmb_read_reference_cache(argv[1],key,2,3,.5,rows,error))return 5;
  }
  cout<<"7 cache round-trip, parameter and integrity checks passed\n";
}
'''
with tempfile.TemporaryDirectory(prefix='qm-cache-') as directory:
    directory=Path(directory);cpp=directory/'test.cpp';exe=directory/'test.exe'
    cpp.write_text(prefix+source[start:stop]+main)
    subprocess.run([args.cxx,'-O2','-std=c++17',str(cpp),'-o',str(exe)],check=True,timeout=120)
    subprocess.run([str(exe),str(directory/'reference.dat')],check=True,timeout=120)
