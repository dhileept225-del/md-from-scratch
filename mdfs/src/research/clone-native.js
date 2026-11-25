// Native LAMMPS (stable 22Jul2025 update 6, gcc) with pair_lj_mine.{h,cpp} cloned from pair_lj_cut; recorded once.
export const NATIVE = {
"hdiff": "--- pair_lj_cut.h\n+++ pair_lj_mine.h\n@@ -13,21 +13,21 @@\n \n #ifdef PAIR_CLASS\n // clang-format off\n-PairStyle(lj/cut,PairLJCut);\n+PairStyle(lj/mine,PairLJMine);\n // clang-format on\n #else\n \n-#ifndef LMP_PAIR_LJ_CUT_H\n-#define LMP_PAIR_LJ_CUT_H\n+#ifndef LMP_PAIR_LJ_MINE_H\n+#define LMP_PAIR_LJ_MINE_H\n \n #include \"pair.h\"\n \n namespace LAMMPS_NS {\n \n-class PairLJCut : public Pair {\n+class PairLJMine : public Pair {\n  public:\n-  PairLJCut(class LAMMPS *);\n-  ~PairLJCut() override;\n+  PairLJMine(class LAMMPS *);\n+  ~PairLJMine() override;\n   void compute(int, int) override;\n   void settings(int, char **) override;\n   void coeff(int, char **) override;\n",
"cppdiff": "--- pair_lj_cut.cpp\n+++ pair_lj_mine.cpp\n@@ -15,7 +15,7 @@\n    Contributing author: Paul Crozier (SNL)\n ------------------------------------------------------------------------- */\n \n-#include \"pair_lj_cut.h\"\n+#include \"pair_lj_mine.h\"\n \n #include \"atom.h\"\n #include \"comm.h\"\n@@ -36,7 +36,7 @@\n \n /* ---------------------------------------------------------------------- */\n \n-PairLJCut::PairLJCut(LAMMPS *lmp) : Pair(lmp)\n+PairLJMine::PairLJMine(LAMMPS *lmp) : Pair(lmp)\n {\n   respa_enable = 1;\n   born_matrix_enable = 1;\n@@ -45,7 +45,7 @@\n \n /* ---------------------------------------------------------------------- */\n \n-PairLJCut::~PairLJCut()\n+PairLJMine::~PairLJMine()\n {\n   if (copymode) return;\n \n@@ -66,7 +66,7 @@\n \n /* ---------------------------------------------------------------------- */\n \n-void PairLJCut::compute(int eflag, int vflag)\n+void PairLJMine::compute(int eflag, int vflag)\n {\n   int i, j, ii, jj, inum, jnum, itype, jtype;\n   double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;\n@@ -140,7 +140,7 @@\n \n /* ---------------------------------------------------------------------- */\n \n-void PairLJCut::compute_inner()\n+void PairLJMine::compute_inner()\n {\n   int i, j, ii, jj, inum, jnum, itype, jtype;\n   double xtmp, ytmp, ztmp, delx, dely, delz, fpair;\n@@ -213,7 +213,7 @@\n \n /* ---------------------------------------------------------------------- */\n \n-void PairLJCut::compute_middle()\n+void PairLJMine::compute_middle()\n {\n   int i, j, ii, jj, inum, jnum, itype, jtype;\n   double xtmp, ytmp, ztmp, delx, dely, delz, fpair;\n@@ -295,7 +295,7 @@\n \n /* ---------------------------------------------------------------------- */\n \n-void PairLJCut::compute_outer(int eflag, int vflag)\n+void PairLJMine::compute_outer(int eflag, int vflag)\n {\n   int i, j, ii, jj, inum, jnum, itype, jtype;\n   double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;\n@@ -394,7 +394,7 @@\n    allocate all arrays\n ------------------------------------------------------------------------- */\n \n-void PairLJCut::allocate()\n+void PairLJMine::allocate()\n {\n   allocated = 1;\n   int n = atom->ntypes + 1;\n@@ -419,7 +419,7 @@\n    global settings\n ------------------------------------------------------------------------- */\n \n-void PairLJCut::settings(int narg, char **arg)\n+void PairLJMine::settings(int narg, char **arg)\n {\n   if (narg != 1) error->all(FLERR, \"Illegal pair_style command\");\n \n@@ -439,7 +439,7 @@\n    set coeffs for one or more type pairs\n ------------------------------------------------------------------------- */\n \n-void PairLJCut::coeff(int narg, char **arg)\n+void PairLJMine::coeff(int narg, char **arg)\n {\n   if (narg < 4 || narg > 5) error->all(FLERR, \"Incorrect args for pair coefficients\" + utils::errorurl(21));\n   if (!allocated) allocate();\n@@ -472,7 +472,7 @@\n    init specific to this pair style\n ------------------------------------------------------------------------- */\n \n-void PairLJCut::init_style()\n+void PairLJMine::init_style()\n {\n   // request regular or rRESPA neighbor list\n \n@@ -498,7 +498,7 @@\n    init for one type pair i,j and corresponding j,i\n ------------------------------------------------------------------------- */\n \n-double PairLJCut::init_one(int i, int j)\n+double PairLJMine::init_one(int i, int j)\n {\n   if (setflag[i][j] == 0) {\n     epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], sigma[i][i], sigma[j][j]);\n@@ -560,7 +560,7 @@\n    proc 0 writes to restart file\n ------------------------------------------------------------------------- */\n \n-void PairLJCut::write_restart(FILE *fp)\n+void PairLJMine::write_restart(FILE *fp)\n {\n   write_restart_settings(fp);\n \n@@ -580,7 +580,7 @@\n    proc 0 reads from restart file, bcasts\n ------------------------------------------------------------------------- */\n \n-void PairLJCut::read_restart(FILE *fp)\n+void PairLJMine::read_restart(FILE *fp)\n {\n   read_restart_settings(fp);\n   allocate();\n@@ -608,7 +608,7 @@\n    proc 0 writes to restart file\n ------------------------------------------------------------------------- */\n \n-void PairLJCut::write_restart_settings(FILE *fp)\n+void PairLJMine::write_restart_settings(FILE *fp)\n {\n   fwrite(&cut_global, sizeof(double), 1, fp);\n   fwrite(&offset_flag, sizeof(int), 1, fp);\n@@ -620,7 +620,7 @@\n    proc 0 reads from restart file, bcasts\n ------------------------------------------------------------------------- */\n \n-void PairLJCut::read_restart_settings(FILE *fp)\n+void PairLJMine::read_restart_settings(FILE *fp)\n {\n   int me = comm->me;\n   if (me == 0) {\n@@ -639,7 +639,7 @@\n    proc 0 writes to data file\n ------------------------------------------------------------------------- */\n \n-void PairLJCut::write_data(FILE *fp)\n+void PairLJMine::write_data(FILE *fp)\n {\n   for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, \"%d %g %g\\n\", i, epsilon[i][i], sigma[i][i]);\n }\n@@ -648,7 +648,7 @@\n    proc 0 writes all pairs to data file\n ------------------------------------------------------------------------- */\n \n-void PairLJCut::write_data_all(FILE *fp)\n+void PairLJMine::write_data_all(FILE *fp)\n {\n   for (int i = 1; i <= atom->ntypes; i++)\n     for (int j = i; j <= atom->ntypes; j++)\n@@ -657,7 +657,7 @@\n \n /* ---------------------------------------------------------------------- */\n \n-double PairLJCut::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,\n+double PairLJMine::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,\n                          double /*factor_coul*/, double factor_lj, double &fforce)\n {\n   double r2inv, r6inv, forcelj, philj;\n@@ -673,7 +673,7 @@\n \n /* ---------------------------------------------------------------------- */\n \n-void PairLJCut::born_matrix(int /*i*/, int /*j*/, int itype, int jtype, double rsq,\n+void PairLJMine::born_matrix(int /*i*/, int /*j*/, int itype, int jtype, double rsq,\n                             double /*factor_coul*/, double factor_lj, double &dupair,\n                             double &du2pair)\n {\n@@ -695,7 +695,7 @@\n \n /* ---------------------------------------------------------------------- */\n \n-void *PairLJCut::extract(const char *str, int &dim)\n+void *PairLJMine::extract(const char *str, int &dim)\n {\n   dim = 2;\n   if (strcmp(str, \"epsilon\") == 0) return (void *) epsilon;\n",
"make": [
"-- GLOB mismatch!",
"-- GLOB mismatch!",
"-- GLOB mismatch!",
"-- Generating style headers...",
"-- Configuring done (2.3s)",
"[  6%] Building CXX object CMakeFiles/lammps.dir/\u2026/src/force.cpp.o",
"[  7%] Building CXX object CMakeFiles/lammps.dir/\u2026/src/lammps.cpp.o",
"[  7%] Building CXX object CMakeFiles/lammps.dir/\u2026/src/pair_lj_mine.cpp.o",
"[  8%] Linking CXX static library liblammps.a",
"[100%] Linking CXX executable lmp"
],
"mine": [
[
"0",
"3.0000000000000018",
"-6.7733680532592473",
"0",
"-2.2744930532592447",
"-3.7033504200856062"
],
[
"50",
"1.6842865267730129",
"-4.8082494418323201",
"0",
"-2.2824512591203403",
"5.5666131491923005"
],
[
"100",
"1.6712577358996483",
"-4.7875608875581506",
"0",
"-2.2813010053596403",
"5.6613913480892046"
],
[
"150",
"1.6444750691309162",
"-4.7471033686005484",
"0",
"-2.2810074430550982",
"5.8614210992427722"
],
[
"200",
"1.6471542413997959",
"-4.7509052858544134",
"0",
"-2.280791606595244",
"5.8805431267470656"
],
[
"250",
"1.6645596668498575",
"-4.7774327356321811",
"0",
"-2.2812174452324636",
"5.7526089005128513"
]
],
"cut": [
[
"0",
"3.0000000000000018",
"-6.7733680532592473",
"0",
"-2.2744930532592447",
"-3.7033504200856062"
],
[
"50",
"1.6842865267730129",
"-4.8082494418323201",
"0",
"-2.2824512591203403",
"5.5666131491923005"
],
[
"100",
"1.6712577358996483",
"-4.7875608875581506",
"0",
"-2.2813010053596403",
"5.6613913480892046"
],
[
"150",
"1.6444750691309162",
"-4.7471033686005484",
"0",
"-2.2810074430550982",
"5.8614210992427722"
],
[
"200",
"1.6471542413997959",
"-4.7509052858544134",
"0",
"-2.280791606595244",
"5.8805431267470656"
],
[
"250",
"1.6645596668498575",
"-4.7774327356321811",
"0",
"-2.2812174452324636",
"5.7526089005128513"
]
],
"loop": {
"cut": "Loop time of 0.3849 on 1 procs for 250 steps with 4000 atoms",
"mine": "Loop time of 0.34066 on 1 procs for 250 steps with 4000 atoms"
},
"table": " for atom types 1 1: i,r,energy,force\n\nLJ\nN 6 R 0.9 2.4\n\n       1  0.9                    6.63611895325291       138.659623994277     \n       2  1.2                   -0.890965287583076     -2.21169334222308     \n       3  1.5                   -0.320336594278575     -1.15802883104616     \n       4  1.8                   -0.11414705558662      -0.368964449086024    \n       5  2.1                   -0.0460946757486002    -0.13014539773546     \n       6  2.4                   -0.0208215955593359    -0.0517801689403695   \n",
"errGuard": "\u2026/src/pair_lj_mine.h:16:19: error: 'PairLJMine' was not declared in this scope"
};
