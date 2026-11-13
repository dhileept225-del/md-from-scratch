import re, os, json, glob
S='/home/claude/lammps-src/lammps/src'
kinds={'PairStyle':'pair_style','BondStyle':'bond_style','AngleStyle':'angle_style','DihedralStyle':'dihedral_style','ImproperStyle':'improper_style','KSpaceStyle':'kspace_style','AtomStyle':'atom_style','FixStyle':'fix','ComputeStyle':'compute','DumpStyle':'dump','RegionStyle':'region','MinimizeStyle':'min_style','IntegrateStyle':'run_style','CommandStyle':'command'}
styles=[]
for h in glob.glob(S+'/**/*.h', recursive=True):
  rel=os.path.relpath(h,S); pkg=rel.split('/')[0] if '/' in rel else ''
  for l in open(h,errors='ignore'):
    m=re.match(r'^(\w+Style)\(\s*([^,]+?)\s*,\s*(.+)\)\s*;?\s*$',l)
    if m and m.group(1) in kinds and not m.group(2).isupper() and 'Deprecated' not in m.group(3):
      cpp=rel[:-2]+'.cpp'
      if not os.path.exists(S+'/'+cpp): cpp=rel
      styles.append([kinds[m.group(1)],m.group(2),m.group(3).strip(),'src/'+cpp,pkg])
styles.sort()
builtin=[]
for i,l in enumerate(open(S+'/input.cpp'),1):
  m=re.search(r'mycmd == "(\w+)"\)',l)
  if m: builtin.append([m.group(1),i])
L=open(S+'/input.cpp').read().split('\n')
for b in builtin:
  fn=re.search(r'mycmd == "\w+"\) (\w+)\(',L[b[1]-1]); fn=fn.group(1) if fn else b[0]
  d=[i for i,l in enumerate(L,1) if l.startswith('void Input::'+fn+'(')]
  b.append(d[0] if d else 0)
sec={'Atom':'atom_style','Integrate':'run_style','Minimize':'min_style','Pair':'pair_style','Bond':'bond_style','Angle':'angle_style','Dihedral':'dihedral_style','Improper':'improper_style','KSpace':'kspace_style','Fix':'fix','Compute':'compute','Region':'region','Dump':'dump'}
br={}; cur=None
for l in open('/tmp/a6/info2.out'):
  m=re.match(r'^(\w+) styles:',l)
  if m: cur=sec.get(m.group(1)); br[cur]=[]; continue
  if l.startswith('Info'): cur=None; continue
  if cur and l.strip(): br[cur]+=l.split()
R='/home/claude/lammps-src/lammps'; top=[]
for d in sorted(os.listdir(R)):
  if d in('build',) or d.startswith('.') or not os.path.isdir(R+'/'+d): continue
  n=sz=0
  for dp,dn,fn in os.walk(R+'/'+d):
    for f in fn: n+=1; sz+=os.lstat(os.path.join(dp,f)).st_size
  top.append([d,n,round(sz/1e6,1)])
core=len(glob.glob(S+'/*.cpp')); pk=sorted(d for d in os.listdir(S) if os.path.isdir(S+'/'+d))
pkn=[[d,len(glob.glob(S+'/'+d+'/*.cpp'))] for d in pk]
json.dump({'tag':'stable_22Jul2025_update6','styles':styles,'builtin':builtin,'browser':br,'top':top,'core':core,'pkgs':pkn},open('/home/claude/md-from-scratch/mdfs/public/data/lammps-src.json','w'),separators=(',',':'))
print(len(styles),len(builtin),{k:len(v) for k,v in br.items()},core,len(pkn)); print(top); print(pkn)
print([s for s in styles if s[1] in('lj/cut','nvt','nve','lj/expand','lj/smooth/linear','create_atoms')])
print(builtin[:12])
