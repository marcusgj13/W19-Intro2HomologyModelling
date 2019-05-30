from pyrosetta import *
from pyrosetta.teaching import *
pyrosetta.init()

p = Pose()
make_pose_from_sequence(p,"A"*10, "fa_standard",auto_termini=True)
p.dump_pdb("unideal_helix.pdb")
score = ScoreFunction()
score.set_weight(fa_atr, 1.0)
score.set_weight(fa_rep, 1.0)
score.set_weight(hbond_lr_bb, 1.0)
score.set_weight(hbond_sr_bb, 1.0)
score.set_weight(hbond_bb_sc, 1.0)
score.set_weight(hbond_sc, 1.0)

for i in range(1, p.total_residue() + 1):
    p.set_phi(i,-60)
    p.set_psi(i,-45.0)

print("Score is " + str(score(p)))

p.dump_pdb("ideal_helix.pdb")
