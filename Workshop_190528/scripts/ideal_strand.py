from pyrosetta import *
from pyrosetta.teaching import *
from rosetta.protocols.relax import *
init()

p = pose_from_sequence("A"*9, "fa_standard")
p.dump_pdb("unideal_sheet.pdb")

score = ScoreFunction()
score.set_weight(fa_atr, 1.0)
score.set_weight(fa_rep, 1.0)
score.set_weight(hbond_lr_bb, 1.0)
score.set_weight(hbond_sr_bb, 1.0)
score.set_weight(hbond_bb_sc, 1.0)
score.set_weight(hbond_sc, 1.0)

#for i in range(1, p.total_residue() + 1):
#    p.set_phi(i, -139)
#    p.set_psi(i, 135)

print("Score is " + str(score(p)))
relax = ClassicRelax()
new_score = get_fa_scorefxn()
relax.set_scorefxn(new_score)
relax.apply(p)
p.dump_pdb("ideal_sheet.pdb")
