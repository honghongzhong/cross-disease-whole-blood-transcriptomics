"""Run the manuscript analyses, then draw the five figures."""
from pathlib import Path
import argparse,os,subprocess,sys,urllib.request
P=Path(__file__).resolve().parent
URL='https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.entrez.gmt'
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--data',type=Path,required=True,help='Directory containing the archived contrasts/ tables');ap.add_argument('--output',type=Path,required=True);ap.add_argument('--hallmark',type=Path,help='Optional already downloaded Hallmark 2025.1.Hs Entrez GMT');args=ap.parse_args()
 out=args.output.resolve()
 if out.exists():raise SystemExit('Choose a new output directory.')
 out.mkdir(parents=True);(out/'resources').mkdir();gmt=out/'resources/h.all.v2025.1.Hs.entrez.gmt'
 if args.hallmark:gmt.write_bytes(args.hallmark.read_bytes())
 else:urllib.request.urlretrieve(URL,gmt)
 if len(gmt.read_text().splitlines())!=50:raise ValueError('Expected the 50-set human Hallmark 2025.1 resource.')
 env=dict(os.environ,BLOOD_DATA_ROOT=str(args.data.resolve()),BLOOD_OUTPUT=str(out))
 for name in ['run_additional_experiments.py','select_dimension.py','finalize_selected_dimensions.py','analyze_residuals.py']:
  subprocess.run([sys.executable,'-u',str(P/'src'/name)],env=env,check=True)
 subprocess.run([sys.executable,str(P/'make_figures.py'),'--results',str(out/'analysis'),'--output',str(out/'figures')],check=True)
 print('Results:',out/'analysis');print('Figures:',out/'figures')
if __name__=='__main__':main()
