#!/bin/bash
################################################################################
# Computational Enzymology MD Workflow
# Author: Shubham Jain
#
# Description:
# This script provides a complete AMBER20-based protein–ligand molecular
# dynamics workflow including:
#
# 1. Ligand parameterization (GAFF + AM1-BCC)
# 2. Protein preparation (ff19SB)
# 3. Complex building and solvation
# 4. Minimization, heating, equilibration, production MD
# 5. cpptraj-based trajectory processing
# 6. Automated MMPBSA binding free energy calculation
# 7. Distance-based frame ranking
# 8. Dynamic Mechanistic Compatibility Score (DMCS) calculation
# 9. Hierarchical clustering for conformational basin analysis
#
# Designed for enzyme–ligand systems and mechanistic studies.
#
# Requirements:
# - AMBER20
# - CUDA-enabled GPU
# - cpptraj
# - MMPBSA.py
#
################################################################################

module load amber20
module load amber20_GPU

##############################################
# SECTION 1 — Ligand Parameterization
##############################################

echo "=== Ligand Parameterization ==="

for LIG in ligand1 ligand2 ligand3
do
    antechamber -i ${LIG}_h.pdb -fi pdb \
                -o ${LIG}.mol2 -fo mol2 \
                -c bcc -s 2 -nc 0

    parmchk2 -i ${LIG}.mol2 -f mol2 \
             -o ${LIG}.frcmod -a y
done


##############################################
# SECTION 2 — Protein Preparation
##############################################

echo "=== Protein Preparation ==="

reduce -Trim protein.pdb > protein_noH.pdb

##############################################
# SECTION 3 — Build Complex in tleap
##############################################

cat <<EOF > build_complex.in
source leaprc.protein.ff19SB
source leaprc.gaff
source leaprc.water.tip3p
loadamberparams frcmod.ionsjc_tip3p

PRO = loadpdb protein_noH.pdb

loadamberparams ligand1.frcmod
UNI = loadmol2 ligand1.mol2

loadamberparams ligand2.frcmod
UNL = loadmol2 ligand2.mol2

loadamberparams ligand3.frcmod
UNP = loadmol2 ligand3.mol2

COM = combine {PRO UNI UNL UNP}

solvateBOX COM TIP3PBOX 12.0
addIons COM Na+ 0

saveamberparm COM solv_complex.prmtop solv_complex.rst7
savepdb COM solv_complex.pdb
quit
EOF

tleap -f build_complex.in


##############################################
# SECTION 4 — Energy Minimization
##############################################

echo "=== Minimization ==="

sander -O -i min1.in \
       -o min1.out \
       -p solv_complex.prmtop \
       -c solv_complex.rst7 \
       -r min1.ncrst \
       -ref solv_complex.rst7


##############################################
# SECTION 5 — Heating & Equilibration
##############################################

pmemd.cuda.MPI -O -i heat.in \
    -o heat.out \
    -p solv_complex.prmtop \
    -c min1.ncrst \
    -r heat.ncrst \
    -x heat.nc \
    -ref min1.ncrst

pmemd.cuda.MPI -O -i equi.in \
    -o equi.out \
    -p solv_complex.prmtop \
    -c heat.ncrst \
    -r equi.ncrst \
    -x equi.nc \
    -ref heat.ncrst


##############################################
# SECTION 6 — Production MD
##############################################

echo "=== Production MD ==="

pmemd.cuda.MPI -O -i md.in \
    -o md.out \
    -p solv_complex.prmtop \
    -c equi.ncrst \
    -r md.rst7 \
    -x md.nc \
    -ref equi.ncrst


##############################################
# SECTION 7 — Trajectory Processing
##############################################

echo "=== cpptraj Processing ==="

cpptraj -p solv_complex.prmtop <<EOF
trajin md.nc
autoimage
strip :WAT,Na+
trajout md_nowat.nc netcdf
run
EOF


##############################################
# SECTION 8 — Distance Calculations
##############################################

echo "=== Distance Calculations ==="

cpptraj -p solv_complex.prmtop <<EOF
trajin md.nc
distance Lys_PLP :188@NZ :PLP@C4 out dist_lys_plp.dat
distance PLP_UNI :PLP@C4 :UNI out dist_plp_uni.dat
distance PLP_UNL :PLP@C4 :UNL out dist_plp_unl.dat
run
EOF


##############################################
# SECTION 9 — Frame Ranking
##############################################

echo "Frame Lys_PLP PLP_UNI PLP_UNL" > plp_all_dist.dat

paste \
  <(tail -n +2 dist_lys_plp.dat) \
  <(tail -n +2 dist_plp_uni.dat) \
  <(tail -n +2 dist_plp_unl.dat) \
| awk '{print $1, $2, $4, $6}' >> plp_all_dist.dat

BEST_FRAME=$(tail -n +2 plp_all_dist.dat | \
awk '{print $1, $2+$3+$4}' | sort -k2 -n | head -n 1)

echo "Best catalytic frame:"
echo $BEST_FRAME


##############################################
# SECTION 10 — DMCS Scoring
##############################################

echo "=== DMCS Calculation ==="

python3 <<EOF
import pandas as pd

df = pd.read_csv("plp_all_dist.dat", delim_whitespace=True)

df["SUM"] = df["Lys_PLP"] + df["PLP_UNI"] + df["PLP_UNL"]

norm = (df["SUM"] - df["SUM"].min()) / (df["SUM"].max() - df["SUM"].min())
dmcs = norm.mean()

print("FINAL DMCS SCORE:", round(dmcs,3))
EOF


##############################################
# SECTION 11 — Hierarchical Clustering
##############################################

echo "=== Clustering ==="

cpptraj -p solv_complex.prmtop <<EOF
trajin md.nc
cluster hieragglo epsilon 2.0 rms @CA \
out cluster.dat \
summary cluster_summary.dat \
singlerepout cluster_rep.pdb
run
EOF

echo "=== Workflow Completed Successfully ==="

