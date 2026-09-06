from pathlib import Path
import runpy,sys
sys.argv[1:1]=['--contrast',Path(__file__).parent.name]
runpy.run_path(str(Path(__file__).parent.parent/'_shared'/'reproduce.py'),run_name='__main__')
