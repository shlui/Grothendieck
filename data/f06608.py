"""Freeze numerical discovery coefficients as exact raw-Hermite rationals.

The result has finite spectral prerequisites only. Complete q and outer-source
certificates are separate; floating discovery values are not proof inputs.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,math,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'lower_amp_certificate'))
from general_R3_cover import validate
from certify_local_point import rational

def freeze(path,b):
 d=json.loads(Path(path).read_text());x=d['history'][-1]['x'];n=d['order'];ia=14+n
 pad=Q(1,10**7)
 while True:
  p=dict(ell=rational(x[0]),beta=rational(x[1]/math.sqrt(6)),p0=rational(d['p'][0]),p2=rational(d['p'][1]/math.sqrt(2)),E00=rational(x[2])+pad,E02=rational(x[3]/math.sqrt(2)),E22=rational(x[4]/2)+pad,R11=rational(x[5])+pad,R13=rational(x[6]/math.sqrt(6)),R33=rational(x[7]/6)+pad,gamma=rational(x[8]),y=rational(x[9]/math.sqrt(2)),r=Q(str(d['r'])),kappa=Q(str(d['kappa'])),eps=Q(b),major=[rational(t)for t in x[10:14]],phi=[rational(t)for t in x[14:14+n]],nuisance0=rational(x[ia]),nuisance2=rational(x[ia+1]),eta=rational(x[ia+2]),zeta=rational(x[ia+3]),c=rational(x[ia+4]),d0=rational(x[ia+5]),G=rational(x[ia+6]))
  try:validate(p);return p,pad
  except AssertionError:
   pad*=10
   if pad>Q(1,100):raise

def main():
 ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('output');ap.add_argument('--b',default='7171/10000');a=ap.parse_args();p,pad=freeze(a.source,a.b)
 Path(a.output).write_text(json.dumps(dict(kind='EXACT_RATIONAL_GENERAL_R3_PROFILE',origin=a.source,rounding='11decimal places, exact positive spectral diagonal repair; no numerical value used for certificate acceptance.',diagonal_repair=str(pad),profile={k:[str(t)for t in v]if isinstance(v,list)else str(v)for k,v in p.items()}),indent=2));print('PASS_EXACT_SIX_SPECTRAL_GATES',str(pad))
if __name__=='__main__':main()
